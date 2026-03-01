"""Pydantic schemas for sequence input and validation."""

from pydantic import BaseModel, Field, field_validator
import re


class SequenceSubmission(BaseModel):
    """Input schema for submitting a biologic sequence for binding prediction."""

    sequence_type: str = Field(
        ...,
        description="Type of biologic: antibody, nanobody, or peptide",
        pattern="^(antibody|nanobody|peptide)$",
    )
    sequence_name: str = Field(
        default="Untitled",
        max_length=128,
        description="User-defined name for this sequence",
    )
    # For nanobodies and peptides — a single chain
    sequence: str | None = Field(
        default=None,
        description="Single amino acid sequence (for nanobody/peptide)",
    )
    # For antibodies — heavy + light chains
    heavy_chain: str | None = Field(
        default=None,
        description="Heavy chain amino acid sequence (for antibody)",
    )
    light_chain: str | None = Field(
        default=None,
        description="Light chain amino acid sequence (for antibody)",
    )

    @field_validator("sequence", "heavy_chain", "light_chain", mode="before")
    @classmethod
    def clean_sequence(cls, v: str | None) -> str | None:
        if v is None:
            return v
        # Remove whitespace and FASTA header lines
        lines = v.strip().split("\n")
        cleaned = "".join(
            line.strip() for line in lines if not line.startswith(">")
        )
        return cleaned.upper()

    @field_validator("sequence", "heavy_chain", "light_chain")
    @classmethod
    def validate_amino_acids(cls, v: str | None) -> str | None:
        if v is None:
            return v
        if not v:
            return v
        valid_pattern = re.compile(r"^[ACDEFGHIKLMNPQRSTVWY]+$")
        if not valid_pattern.match(v):
            invalid_chars = set(v) - set("ACDEFGHIKLMNPQRSTVWY")
            raise ValueError(
                f"Invalid amino acid characters: {invalid_chars}. "
                "Only standard 20 amino acids are accepted."
            )
        if len(v) < 10:
            raise ValueError("Sequence must be at least 10 amino acids long.")
        if len(v) > 5000:
            raise ValueError("Sequence must be at most 5000 amino acids long.")
        return v

    def get_primary_sequence(self) -> str:
        """Get the primary sequence for this submission."""
        if self.sequence_type == "antibody":
            parts = []
            if self.heavy_chain:
                parts.append(self.heavy_chain)
            if self.light_chain:
                parts.append(self.light_chain)
            return ":".join(parts) if parts else ""
        return self.sequence or ""


class SequenceValidationResult(BaseModel):
    """Result of sequence validation."""

    valid: bool
    sequence_length: int = 0
    sequence_type: str = ""
    warnings: list[str] = []
    errors: list[str] = []
