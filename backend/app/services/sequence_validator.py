"""Sequence validation utilities for biologic inputs."""

import re
from app.schemas.sequence import SequenceSubmission, SequenceValidationResult


def validate_submission(submission: SequenceSubmission) -> SequenceValidationResult:
    """Validate a sequence submission and return detailed results."""
    errors: list[str] = []
    warnings: list[str] = []
    seq_length = 0

    if submission.sequence_type == "antibody":
        if not submission.heavy_chain and not submission.light_chain:
            errors.append("Antibody submissions require at least a heavy chain sequence.")
        if submission.heavy_chain:
            seq_length += len(submission.heavy_chain)
            if len(submission.heavy_chain) < 100:
                warnings.append(
                    f"Heavy chain is unusually short ({len(submission.heavy_chain)} aa). "
                    "Typical heavy chains are 440-470 aa."
                )
            if len(submission.heavy_chain) > 600:
                warnings.append(
                    f"Heavy chain is unusually long ({len(submission.heavy_chain)} aa)."
                )
        if submission.light_chain:
            seq_length += len(submission.light_chain)
            if len(submission.light_chain) < 90:
                warnings.append(
                    f"Light chain is unusually short ({len(submission.light_chain)} aa). "
                    "Typical light chains are 210-230 aa."
                )
        if not submission.light_chain:
            warnings.append(
                "No light chain provided. Predictions may be less accurate "
                "without the complete antibody structure."
            )
    elif submission.sequence_type == "nanobody":
        if not submission.sequence:
            errors.append("Nanobody submissions require a sequence.")
        else:
            seq_length = len(submission.sequence)
            if seq_length < 100:
                warnings.append(
                    f"Nanobody sequence is short ({seq_length} aa). "
                    "Typical VHH domains are 110-130 aa."
                )
            if seq_length > 200:
                warnings.append(
                    f"Nanobody sequence is long ({seq_length} aa). "
                    "Consider if this is a single VHH domain."
                )
    elif submission.sequence_type == "peptide":
        if not submission.sequence:
            errors.append("Peptide submissions require a sequence.")
        else:
            seq_length = len(submission.sequence)
            if seq_length > 100:
                warnings.append(
                    f"Peptide sequence is long ({seq_length} aa). "
                    "Consider classifying as protein instead."
                )

    return SequenceValidationResult(
        valid=len(errors) == 0,
        sequence_length=seq_length,
        sequence_type=submission.sequence_type,
        warnings=warnings,
        errors=errors,
    )
