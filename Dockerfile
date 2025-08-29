# ---- Build Stage ----
FROM python:3.13-slim AS builder
COPY --from=ghcr.io/astral-sh/uv:latest /uv /uvx /bin/

WORKDIR /app

RUN python -m venv /opt/venv

COPY tropic-core ./tropic-core/
COPY tropic-api ./tropic-api/

ENV PATH="/opt/venv/bin:$PATH"
RUN uv pip install ./tropic-core ./tropic-api[strict]

# ---- Run Stage ----
FROM python:3.13-slim

WORKDIR /app

COPY --from=builder /opt/venv /opt/venv
COPY data/db.pkl ./db.pkl

ENV \
  PYTHONDONTWRITEBYTECODE=1 \
  PYTHONUNBUFFERED=1 \
  PATH="/opt/venv/bin:${PATH}"

EXPOSE 8000
CMD ["/bin/sh", "-c", "tropic-load --filename db.pkl && tropic-api"]
