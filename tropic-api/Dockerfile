# Dockerfile for a FastAPI application with a multi-stage build
# This Dockerfile uses a multi-stage build to create a lightweight image for a FastAPI application.

# Stage 1: Builder
FROM python:3.13-slim as builder

WORKDIR /app

# Install build dependencies
RUN apt-get update && \
    apt-get install -y --no-install-recommends gcc python3-dev

# Create and activate virtual environment
RUN python -m venv /opt/venv
ENV PATH="/opt/venv/bin:$PATH"

# Install dependencies
COPY requirements.txt .
RUN pip install --no-cache-dir --upgrade pip && \
    pip install --no-cache-dir .[strict]

# Stage 2: Runtime
FROM python:3.13-slim

WORKDIR /app

# Copy virtual environment from builder
COPY --from=builder /opt/venv /opt/venv
ENV PATH="/opt/venv/bin:$PATH"

# Install runtime dependencies
RUN apt-get update && \
    apt-get install -y --no-install-recommends libgomp1 && \
    apt-get clean && \
    rm -rf /var/lib/apt/lists/*

# Copy application files
COPY main.py .

EXPOSE 8000
CMD ["tropic-api"]
