# HumanProof — Deployment Guide

Primary deployment targets: **Vercel** (frontend) + **Railway** or **Google Cloud Run** (backend).
Docker/VPS and GPU worker sections are retained at the bottom for future reference.

---

## Table of Contents

1. [Local Development](#1-local-development)
2. [Environment Variables Reference](#2-environment-variables-reference)
3. [Frontend → Vercel](#3-frontend--vercel)
4. [Backend → Railway](#4-backend--railway-recommended)
5. [Backend → Google Cloud Run](#5-backend--google-cloud-run)
6. [PostgreSQL Migration](#6-postgresql-migration)
7. [Loading Real Data in Production](#7-loading-real-data-in-production)
8. [Health Checks & Monitoring](#8-health-checks--monitoring)
9. [Troubleshooting](#9-troubleshooting)
10. [Future: Docker Deployment (Self-Hosted VPS)](#10-future-docker-deployment-self-hosted-vps)
11. [Future: Nginx Reverse Proxy & SSL](#11-future-nginx-reverse-proxy--ssl)
12. [Future: GPU Worker Setup (Boltz-2 / Phase 3)](#12-future-gpu-worker-setup-boltz-2--phase-3)

---

## 1. Local Development

### Prerequisites

| Tool | Version | Install |
|------|---------|---------|
| Python | ≥ 3.11 | [python.org](https://python.org) or `conda` |
| Node.js | ≥ 20 LTS | [nodejs.org](https://nodejs.org) |
| npm | ≥ 10 | bundled with Node.js |

### Backend

```bash
cd backend

# Create a virtual environment
python3 -m venv .venv
source .venv/bin/activate        # Windows: .venv\Scripts\activate

# Install dependencies
pip install -e ".[dev]"

# Populate the database (run once; re-run to reset)
python load_real_data.py

# Start the API server
uvicorn app.main:app --reload --host 127.0.0.1 --port 8000
```

The API is now available at `http://127.0.0.1:8000`.
Interactive docs: `http://127.0.0.1:8000/docs`

### Frontend

```bash
cd frontend

# Install dependencies
npm install

# Start the dev server (proxies /api/* → localhost:8000 automatically)
npm run dev
```

The app is now available at `http://localhost:3000`.

### Verify the stack

```bash
curl http://127.0.0.1:8000/health
# → {"status":"ok","app":"HumanProof"}

curl "http://127.0.0.1:8000/api/v1/targets/search?q=TP5"
# → [{"gene_symbol":"TP53", ...}]
```

---

## 2. Environment Variables Reference

### Backend — `HUMANPROOF_*` prefix (pydantic-settings)

| Variable | Default | Description |
|----------|---------|-------------|
| `HUMANPROOF_APP_NAME` | `HumanProof` | App display name |
| `HUMANPROOF_DATABASE_URL` | `sqlite+aiosqlite:///./humanproof.db` | SQLAlchemy async database URL |
| `HUMANPROOF_CORS_ORIGINS` | `["http://localhost:3000"]` | JSON list of allowed CORS origins |
| `HUMANPROOF_DATA_DIR` | `<repo>/data` | Path to data directory (SHAP JSON, pickles) |

Set via environment or a `.env` file in `backend/`:

```bash
# backend/.env  (never commit this file)
HUMANPROOF_DATABASE_URL=postgresql+asyncpg://user:pass@db:5432/humanproof
HUMANPROOF_CORS_ORIGINS=["https://your-project.vercel.app"]
HUMANPROOF_DATA_DIR=/app/data
```

### Frontend — `NEXT_PUBLIC_*` prefix

| Variable | Default | Description |
|----------|---------|-------------|
| `NEXT_PUBLIC_API_URL` | `http://localhost:8000` | Backend base URL (used by browser in prod) |

> **Note:** In development the Next.js rewrite rule in `next.config.ts` proxies all `/api/*`
> requests to `http://127.0.0.1:8000`, so `NEXT_PUBLIC_API_URL` is only needed in production
> when the frontend and backend are on separate domains.

---

## 3. Frontend → Vercel

Vercel detects Next.js automatically and handles builds, previews, and CDN globally.

### 3.1 Deploy

1. Go to [vercel.com/new](https://vercel.com/new) and import the GitHub repo.
2. Set **Root Directory** to `frontend/`.
3. Add the environment variable:
   ```
   NEXT_PUBLIC_API_URL=https://<your-backend-url>
   ```
   (Use the Railway or Cloud Run URL from sections 4 or 5.)
4. Click **Deploy**. Vercel auto-detects Next.js; no build command changes needed.

### 3.2 Custom domain (optional)

In the Vercel dashboard → **Domains** → add your domain. Vercel provisions SSL automatically.

### 3.3 Update CORS after deploy

Once you have your Vercel production URL (e.g., `https://humanproof.vercel.app`), update
`HUMANPROOF_CORS_ORIGINS` on the backend to include it:

```bash
HUMANPROOF_CORS_ORIGINS=["https://humanproof.vercel.app"]
```

### 3.4 Redeployments

Push to `main` triggers an automatic redeploy. Feature branches get preview URLs automatically.

---

## 4. Backend → Railway (Recommended)

Railway is the simplest managed backend: free tier, automatic deploys from GitHub,
built-in PostgreSQL, persistent volumes for large data files.

### 4.1 Create project

1. Go to [railway.app](https://railway.app) → **New Project** → **Deploy from GitHub repo**.
2. Select the repo and set **Root Directory** to `backend/`.

### 4.2 Add PostgreSQL

In the Railway project dashboard → **Add Service** → **Database** → **PostgreSQL**.
Railway auto-injects `DATABASE_URL` (postgres:// format) into all services.

### 4.3 Procfile

Create `backend/Procfile`:

```
web: uvicorn app.main:app --host 0.0.0.0 --port $PORT --workers 1
```

Use `--workers 1` when running with SQLite (single-file DB) to avoid startup race conditions.

### 4.4 Environment variables

In the Railway backend service → **Variables**, add:

```
HUMANPROOF_DATABASE_URL=${{Postgres.DATABASE_URL}}
HUMANPROOF_CORS_ORIGINS=["https://humanproof.vercel.app"]
HUMANPROOF_DATA_DIR=/app/data
```

> **asyncpg URL:** Railway's `DATABASE_URL` uses `postgres://` scheme. FastAPI needs
> `postgresql+asyncpg://`. Add a startup shim or set a separate variable:
>
> ```bash
> # In Railway Variables panel, add:
> HUMANPROOF_DATABASE_URL=postgresql+asyncpg://${{Postgres.PGUSER}}:${{Postgres.PGPASSWORD}}@${{Postgres.PGHOST}}:${{Postgres.PGPORT}}/${{Postgres.PGDATABASE}}
> ```

### 4.5 Persistent volume for data files

The SHAP JSON (`gene_shap_dr.json`, 455 MB), pLoF pickle, and CellxGene CSVs must be
accessible at runtime but are too large to commit to git.

1. In Railway → **Add Volume** → mount at `/app/data`.
2. Sync files from Zenodo (recommended):

```bash
# Install Railway CLI
npm install -g @railway/cli
railway login
railway link   # link to your project

# Canonical sync workflow (run from repo root)
# Defaults to Zenodo record 18827087 and verifies expected file sizes
./sync_data_from_urls.sh
```

Optional: pin a specific dataset record/version:

```bash
ZENODO_RECORD_ID=18827087 ./sync_data_from_urls.sh
```

> Note: `railway ssh` stdin streaming and `railway volume cp` may be unreliable or unavailable
> depending on Railway CLI version. `sync_data_from_urls.sh` is the supported method.

### 4.6 Seed the database

After the first deploy, seed from inside the running service:

```bash
railway ssh "cd /app && /app/.venv/bin/python load_real_data.py"
```

This populates PostgreSQL from the uploaded data files (~19K genes, takes several minutes).

### 4.7 Get the backend URL

Railway dashboard → backend service → **Settings** → **Networking** → copy the public URL
(e.g., `https://humanproof-backend.up.railway.app`). Use this as `NEXT_PUBLIC_API_URL` in Vercel.

### 4.8 Redeployments

Push to `main` triggers an automatic redeploy. The volume persists across deploys.

---

## 5. Backend → Google Cloud Run

Cloud Run is a good alternative if you want Google's infrastructure, pay-per-request billing,
or need to stay within a Google Cloud organization.

### 5.1 Prerequisites

```bash
# Install gcloud CLI: https://cloud.google.com/sdk/docs/install
gcloud auth login
gcloud config set project YOUR_PROJECT_ID
gcloud services enable run.googleapis.com sqladmin.googleapis.com artifactregistry.googleapis.com
```

### 5.2 Create Cloud SQL (PostgreSQL)

```bash
gcloud sql instances create humanproof-db \
  --database-version=POSTGRES_16 \
  --tier=db-f1-micro \
  --region=us-central1

gcloud sql databases create humanproof --instance=humanproof-db

gcloud sql users set-password postgres \
  --instance=humanproof-db \
  --password=CHANGE_ME
```

Note your instance connection name (format: `PROJECT:REGION:INSTANCE`):
```bash
gcloud sql instances describe humanproof-db --format="value(connectionName)"
```

### 5.3 Store data files in Cloud Storage

```bash
# Create a bucket
gsutil mb gs://humanproof-data

# Upload large data files
gsutil cp data/safety_model_output/dr/gene_shap_dr.json gs://humanproof-data/safety_model_output/dr/
gsutil cp data/genebass_pLoF_filtered.pkl gs://humanproof-data/
gsutil cp data/LOEUF_scores.csv.gz gs://humanproof-data/
```

> **Option A (simpler):** Download from Cloud Storage on container startup using a startup script.
> **Option B:** Include smaller data files directly in the Docker image; download the 455 MB SHAP JSON at startup.

Startup script `backend/start.sh`:

```bash
#!/bin/bash
set -e

# Download data files from GCS if not already present
mkdir -p /app/data/safety_model_output/dr
gsutil -q cp gs://humanproof-data/safety_model_output/dr/gene_shap_dr.json \
  /app/data/safety_model_output/dr/gene_shap_dr.json 2>/dev/null || true
gsutil -q cp gs://humanproof-data/genebass_pLoF_filtered.pkl /app/data/ 2>/dev/null || true
gsutil -q cp gs://humanproof-data/LOEUF_scores.csv.gz /app/data/ 2>/dev/null || true

exec uvicorn app.main:app --host 0.0.0.0 --port "$PORT" --workers 2
```

### 5.4 Dockerfile for Cloud Run

Create `docker/Dockerfile.cloudrun`:

```dockerfile
FROM python:3.11-slim

WORKDIR /app

RUN apt-get update && apt-get install -y --no-install-recommends \
    gcc libpq-dev curl \
    && rm -rf /var/lib/apt/lists/*

# Install gcloud storage component (lightweight gsutil alternative)
RUN pip install --no-cache-dir google-cloud-storage

COPY backend/pyproject.toml ./
RUN pip install --no-cache-dir ".[standard]" "asyncpg>=0.29.0"

COPY backend/app ./app
COPY docker/start.sh ./start.sh
RUN chmod +x start.sh

EXPOSE 8080

CMD ["./start.sh"]
```

### 5.5 Build and push the image

```bash
# Configure Docker to use Artifact Registry
gcloud auth configure-docker us-central1-docker.pkg.dev

# Create repository
gcloud artifacts repositories create humanproof \
  --repository-format=docker \
  --location=us-central1

# Build and push
docker build -f docker/Dockerfile.cloudrun \
  -t us-central1-docker.pkg.dev/YOUR_PROJECT_ID/humanproof/backend:latest .
docker push us-central1-docker.pkg.dev/YOUR_PROJECT_ID/humanproof/backend:latest
```

### 5.6 Deploy to Cloud Run

```bash
gcloud run deploy humanproof-backend \
  --image us-central1-docker.pkg.dev/YOUR_PROJECT_ID/humanproof/backend:latest \
  --region us-central1 \
  --platform managed \
  --allow-unauthenticated \
  --memory 2Gi \
  --cpu 2 \
  --timeout 300 \
  --add-cloudsql-instances YOUR_PROJECT_ID:us-central1:humanproof-db \
  --set-env-vars "HUMANPROOF_DATABASE_URL=postgresql+asyncpg://postgres:CHANGE_ME@/humanproof?host=/cloudsql/YOUR_PROJECT_ID:us-central1:humanproof-db" \
  --set-env-vars "HUMANPROOF_CORS_ORIGINS=[\"https://humanproof.vercel.app\"]" \
  --set-env-vars "HUMANPROOF_DATA_DIR=/app/data"
```

### 5.7 Seed the database

Cloud Run does not have an interactive shell by default. Run seeding as a one-off job:

```bash
gcloud run jobs create humanproof-seed \
  --image us-central1-docker.pkg.dev/YOUR_PROJECT_ID/humanproof/backend:latest \
  --region us-central1 \
  --add-cloudsql-instances YOUR_PROJECT_ID:us-central1:humanproof-db \
  --set-env-vars "HUMANPROOF_DATABASE_URL=postgresql+asyncpg://postgres:CHANGE_ME@/humanproof?host=/cloudsql/YOUR_PROJECT_ID:us-central1:humanproof-db" \
  --set-env-vars "HUMANPROOF_DATA_DIR=/app/data" \
  --command "python,load_real_data.py"

gcloud run jobs execute humanproof-seed --region us-central1 --wait
```

### 5.8 Get the backend URL

```bash
gcloud run services describe humanproof-backend \
  --region us-central1 \
  --format="value(status.url)"
```

Use this URL as `NEXT_PUBLIC_API_URL` in Vercel and update `HUMANPROOF_CORS_ORIGINS` accordingly.

---

## 6. PostgreSQL Migration

SQLite is used in development. For Railway or Cloud Run deployments, switch to PostgreSQL.

### 6.1 Install the async PostgreSQL driver

```bash
# In backend/ venv
pip install asyncpg psycopg2-binary
```

Add to `pyproject.toml` optional dependencies:

```toml
[project.optional-dependencies]
postgres = ["asyncpg>=0.29.0", "psycopg2-binary>=2.9.0"]
```

### 6.2 Set the database URL

```bash
# PostgreSQL async URL format
HUMANPROOF_DATABASE_URL=postgresql+asyncpg://user:password@host:5432/humanproof
```

The rest of the code (`database.py`, all models) is already compatible — SQLAlchemy 2.0 handles both SQLite and PostgreSQL identically.

### 6.3 Migrate existing SQLite data (if needed)

```bash
# Install pgloader
brew install pgloader     # macOS
apt install pgloader      # Ubuntu

# Run migration
pgloader sqlite:///backend/humanproof.db \
  postgresql://user:password@host:5432/humanproof
```

### 6.4 Production DB sizing guide

| Scale | Setup | Cost |
|-------|-------|------|
| Lab / demo (≤10 users) | Railway Postgres (1 GB RAM) | ~$5/mo |
| Department (≤100 users) | Railway Postgres (2 GB RAM) | ~$25/mo |
| Multi-lab (>100 users) | Cloud SQL db.g1-small or RDS t3.medium | ~$30–50/mo |

---

## 7. Loading Real Data in Production

The backend reads from three large data sources at startup or on first query:
- `gene_shap_dr.json` (455 MB) — per-gene SHAP values, read by `targets.py`
- `genebass_pLoF_filtered.pkl` — pLoF associations, used during DB seeding
- `LOEUF_scores.csv.gz` — constraint scores, used during DB seeding

### 7.1 Railway: upload via CLI or volume

```bash
# Link your local repo to Railway project
railway login && railway link

# Sync all required files to mounted volume (/app/data)
# Uses Zenodo API content URLs and verifies expected byte sizes
./sync_data_from_urls.sh

# Run DB seeding
railway ssh "cd /app && /app/.venv/bin/python load_real_data.py"
```

### 7.2 Google Cloud Run: upload to Cloud Storage

```bash
# Upload / update data files
gsutil cp data/safety_model_output/dr/gene_shap_dr.json \
  gs://humanproof-data/safety_model_output/dr/
gsutil cp data/genebass_pLoF_filtered.pkl gs://humanproof-data/
gsutil cp data/LOEUF_scores.csv.gz gs://humanproof-data/

# Re-run the seed job (drops and recreates DB tables)
gcloud run jobs execute humanproof-seed --region us-central1 --wait
```

### 7.3 Data update workflow (both platforms)

```bash
# 1. Upload updated data files (see 7.1 or 7.2)

# 2. Backup existing DB first (Railway example)
railway run pg_dump $DATABASE_URL > backup_$(date +%Y%m%d).sql

# 3. Re-seed
railway ssh "cd /app && /app/.venv/bin/python load_real_data.py"

# 4. Verify
curl https://<backend-url>/api/v1/targets/search?q=BRCA1
```

---

## 8. Health Checks & Monitoring

### 8.1 Built-in health endpoint

```bash
curl https://<backend-url>/health
# → {"status": "ok", "app": "HumanProof"}
```

### 8.2 Uptime monitoring

Use [UptimeRobot](https://uptimerobot.com) (free) or [Betterstack](https://betterstack.com):

- Monitor URL: `https://<backend-url>/health`
- Check interval: 5 minutes
- Alert channels: email / Slack

### 8.3 Railway logs

```bash
railway logs          # tail live logs
railway logs --tail 200  # last 200 lines
```

### 8.4 Google Cloud Run logs

```bash
gcloud run services logs read humanproof-backend \
  --region us-central1 --limit 100
```

Or view in **Cloud Console** → Cloud Run → humanproof-backend → Logs.

---

## 9. Troubleshooting

### Backend won't start: "database is locked"

SQLite doesn't support concurrent writes. Either:
- Run a single uvicorn worker (`--workers 1`), or
- Switch to PostgreSQL (required for Railway/Cloud Run multi-instance deployments)

### CORS errors in browser

Add the production frontend URL to `HUMANPROOF_CORS_ORIGINS`:
```bash
HUMANPROOF_CORS_ORIGINS=["https://humanproof.vercel.app","https://www.humanproof.vercel.app"]
```

Note: the value must be a valid JSON array string.

### Frontend shows "Failed to fetch" / blank data

1. Confirm the backend is running: `curl https://<backend-url>/health`
2. Confirm `NEXT_PUBLIC_API_URL` is set correctly in the Vercel environment variables
3. Check browser console — the URL may still point to `localhost:8000` if the env var was missing at build time (trigger a redeploy after adding it)

### asyncpg URL format error (Railway)

Railway's `${{Postgres.DATABASE_URL}}` uses `postgres://` scheme; FastAPI requires `postgresql+asyncpg://`. Set the variable manually using the individual Postgres credential references:

```
HUMANPROOF_DATABASE_URL=postgresql+asyncpg://${{Postgres.PGUSER}}:${{Postgres.PGPASSWORD}}@${{Postgres.PGHOST}}:${{Postgres.PGPORT}}/${{Postgres.PGDATABASE}}
```

### SHAP JSON not loading (gene risk scores missing)

The `gene_shap_dr.json` (455 MB) must be present at `$HUMANPROOF_DATA_DIR/safety_model_output/dr/gene_shap_dr.json`.
Confirm the file is uploaded and the env var points to the right directory.

### Database tables missing after deploy

The app creates tables automatically via `init_db()`. If tables are missing, trigger re-initialization:

```bash
# Railway
railway ssh "cd /app && /app/.venv/bin/python load_real_data.py"

# Cloud Run
gcloud run jobs execute humanproof-seed --region us-central1 --wait
```

### Port 8000 already in use (local dev)

```bash
lsof -ti:8000 | xargs kill -9
```

### pLoF pickle file not found

The generator looks for the file at `$HUMANPROOF_DATA_DIR/genebass_pLoF_filtered.pkl`.
Set the env var or ensure the file is uploaded to the data volume/bucket.

---

## 10. Future: Docker Deployment (Self-Hosted VPS)

> **When to use:** If you need full infrastructure control, a private network, or want to
> self-host on a university HPC or cloud VM.

### Dockerfile — Backend

Create `docker/Dockerfile.backend`:

```dockerfile
FROM python:3.11-slim

WORKDIR /app

RUN apt-get update && apt-get install -y --no-install-recommends \
    gcc libpq-dev curl \
    && rm -rf /var/lib/apt/lists/*

COPY backend/pyproject.toml ./
RUN pip install --no-cache-dir ".[standard]" \
    "asyncpg>=0.29.0" \
    "psycopg2-binary>=2.9.0"

COPY backend/app ./app
COPY data ./data

HEALTHCHECK --interval=30s --timeout=10s --start-period=15s --retries=3 \
  CMD curl -f http://localhost:8000/health || exit 1

EXPOSE 8000

CMD ["uvicorn", "app.main:app", "--host", "0.0.0.0", "--port", "8000", \
     "--workers", "2", "--proxy-headers", "--forwarded-allow-ips=*"]
```

### Dockerfile — Frontend

Create `docker/Dockerfile.frontend`:

```dockerfile
FROM node:20-alpine AS deps
WORKDIR /app
COPY frontend/package*.json ./
RUN npm ci --omit=dev

FROM node:20-alpine AS builder
WORKDIR /app
COPY --from=deps /app/node_modules ./node_modules
COPY frontend/ .

ARG NEXT_PUBLIC_API_URL=https://api.humanproof.yourdomain.com
ENV NEXT_PUBLIC_API_URL=$NEXT_PUBLIC_API_URL

RUN npm run build

FROM node:20-alpine AS runner
WORKDIR /app
ENV NODE_ENV=production
COPY --from=builder /app/.next/standalone ./
COPY --from=builder /app/.next/static ./.next/static
COPY --from=builder /app/public ./public

EXPOSE 3000
CMD ["node", "server.js"]
```

> Add `output: "standalone"` to `frontend/next.config.ts` to enable the optimized build.

### Docker Compose

Create `docker/docker-compose.yml`:

```yaml
version: "3.9"

services:
  db:
    image: postgres:16-alpine
    restart: unless-stopped
    environment:
      POSTGRES_USER: humanproof
      POSTGRES_PASSWORD: ${DB_PASSWORD}
      POSTGRES_DB: humanproof
    volumes:
      - pg_data:/var/lib/postgresql/data
    healthcheck:
      test: ["CMD-SHELL", "pg_isready -U humanproof"]
      interval: 10s
      timeout: 5s
      retries: 5

  backend:
    build:
      context: ..
      dockerfile: docker/Dockerfile.backend
    restart: unless-stopped
    depends_on:
      db:
        condition: service_healthy
    environment:
      HUMANPROOF_DATABASE_URL: postgresql+asyncpg://humanproof:${DB_PASSWORD}@db:5432/humanproof
      HUMANPROOF_CORS_ORIGINS: '["https://${DOMAIN}"]'
      HUMANPROOF_DATA_DIR: /app/data
    volumes:
      - ../data:/app/data:ro
    ports:
      - "8000:8000"

  frontend:
    build:
      context: ..
      dockerfile: docker/Dockerfile.frontend
      args:
        NEXT_PUBLIC_API_URL: https://api.${DOMAIN}
    restart: unless-stopped
    depends_on:
      - backend
    ports:
      - "3000:3000"

  nginx:
    image: nginx:alpine
    restart: unless-stopped
    depends_on:
      - frontend
      - backend
    ports:
      - "80:80"
      - "443:443"
    volumes:
      - ./nginx.conf:/etc/nginx/conf.d/default.conf:ro
      - certbot_www:/var/www/certbot:ro
      - certbot_certs:/etc/letsencrypt:ro

  certbot:
    image: certbot/certbot
    volumes:
      - certbot_www:/var/www/certbot
      - certbot_certs:/etc/letsencrypt
    entrypoint: >
      sh -c "trap exit TERM; while :; do certbot renew --webroot
      -w /var/www/certbot --quiet; sleep 12h & wait $${!}; done"

volumes:
  pg_data:
  certbot_www:
  certbot_certs:
```

Create `docker/.env` (never commit):

```bash
DOMAIN=humanproof.yourdomain.com
DB_PASSWORD=change_me_to_a_strong_password
```

### Build and launch

```bash
cd docker
docker compose up -d --build
docker compose exec backend python load_real_data.py
docker compose logs -f backend frontend
```

---

## 11. Future: Nginx Reverse Proxy & SSL

> **When to use:** With the Docker/VPS deployment above. Vercel and Cloud Run handle SSL automatically; Railway handles SSL via its proxy.

Create `docker/nginx.conf`:

```nginx
server {
    listen 80;
    server_name humanproof.yourdomain.com api.humanproof.yourdomain.com;

    location /.well-known/acme-challenge/ {
        root /var/www/certbot;
    }

    location / {
        return 301 https://$host$request_uri;
    }
}

server {
    listen 443 ssl http2;
    server_name humanproof.yourdomain.com;

    ssl_certificate     /etc/letsencrypt/live/humanproof.yourdomain.com/fullchain.pem;
    ssl_certificate_key /etc/letsencrypt/live/humanproof.yourdomain.com/privkey.pem;
    ssl_protocols       TLSv1.2 TLSv1.3;

    add_header Strict-Transport-Security "max-age=63072000" always;
    add_header X-Frame-Options DENY;
    add_header X-Content-Type-Options nosniff;

    gzip on;
    gzip_types text/plain text/css application/json application/javascript;

    location / {
        proxy_pass         http://frontend:3000;
        proxy_http_version 1.1;
        proxy_set_header   Upgrade $http_upgrade;
        proxy_set_header   Connection 'upgrade';
        proxy_set_header   Host $host;
        proxy_cache_bypass $http_upgrade;
    }
}

server {
    listen 443 ssl http2;
    server_name api.humanproof.yourdomain.com;

    ssl_certificate     /etc/letsencrypt/live/humanproof.yourdomain.com/fullchain.pem;
    ssl_certificate_key /etc/letsencrypt/live/humanproof.yourdomain.com/privkey.pem;
    ssl_protocols       TLSv1.2 TLSv1.3;

    add_header Strict-Transport-Security "max-age=63072000" always;

    client_max_body_size 10M;

    limit_req_zone $binary_remote_addr zone=api:10m rate=30r/m;
    limit_req zone=api burst=10 nodelay;

    location / {
        proxy_pass         http://backend:8000;
        proxy_http_version 1.1;
        proxy_set_header   Host $host;
        proxy_set_header   X-Real-IP $remote_addr;
        proxy_set_header   X-Forwarded-For $proxy_add_x_forwarded_for;
        proxy_set_header   X-Forwarded-Proto $scheme;
        proxy_read_timeout 120s;
    }
}
```

---

## 12. Future: GPU Worker Setup (Boltz-2 / Phase 3)

> **When to use:** When swapping `MockBindingPredictor` for real structure-based predictions
> (Boltz-2, AlphaFold-Multimer, ESMFold).

### Architecture

```
Frontend → FastAPI (CPU, Railway or Cloud Run)
               ↓ submits job
         Redis (job queue)
               ↓ picks up job
         Celery Worker (GPU instance)
               ↓ writes results
         PostgreSQL
```

### Add Celery + Redis

```bash
pip install celery redis
```

### GPU instance recommendations

| Use case | Instance | VRAM | Cost |
|----------|----------|------|------|
| Development/testing | NVIDIA RTX 4090 (local) | 24 GB | hardware |
| Small lab (on-demand) | AWS g4dn.xlarge (T4) | 16 GB | ~$0.53/hr |
| Production throughput | AWS g5.xlarge (A10G) | 24 GB | ~$1.01/hr |

### Swap in real predictor

Edit `backend/app/services/binding_service.py`:

```python
def get_predictor() -> BindingPredictor:
    predictor_type = os.getenv("BINDING_PREDICTOR", "mock")
    if predictor_type == "boltz2":
        from app.services.boltz2_predictor import Boltz2Predictor
        return Boltz2Predictor()
    elif predictor_type == "esm":
        from app.services.esm_predictor import ESMFoldPredictor
        return ESMFoldPredictor()
    return MockBindingPredictor()
```

Set `BINDING_PREDICTOR=boltz2` on the GPU worker.

---

## Quick Reference

```
Local dev:
  Backend  →  cd backend && uvicorn app.main:app --reload --port 8000
  Frontend →  cd frontend && npm run dev
  Seed DB  →  cd backend && python load_real_data.py

Vercel + Railway:
  Frontend →  vercel.com/new → root: frontend/ → set NEXT_PUBLIC_API_URL
  Backend  →  railway.app → root: backend/ → add Postgres plugin → set env vars
  Data     →  ./sync_data_from_urls.sh
  Seed     →  railway ssh "cd /app && /app/.venv/bin/python load_real_data.py"

Vercel + Google Cloud Run:
  Frontend →  vercel.com/new → root: frontend/ → set NEXT_PUBLIC_API_URL
  Backend  →  gcloud run deploy humanproof-backend --image ...
  Data     →  gsutil cp data/ gs://humanproof-data/
  Seed     →  gcloud run jobs execute humanproof-seed --wait

Key URLs:
  App      →  https://humanproof.vercel.app
  API      →  https://<railway-or-cloudrun-url>
  API docs →  https://<railway-or-cloudrun-url>/docs
  Health   →  https://<railway-or-cloudrun-url>/health
```
