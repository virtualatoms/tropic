# Stage 1: Builder
FROM python:3.13-slim

WORKDIR /app

RUN pip install uv

# Copy the entire project context (both tropic-core and tropic-api) into the container
# This is necessary so uv can find and install both local packages
COPY tropic-core ./tropic-core/
COPY tropic-api ./tropic-api/

# Use uv to install both local packages in editable mode.
# The `-e` flag is important for monorepos as it correctly links the local packages.
# This command will read the pyproject.toml from each directory and install dependencies.
RUN uv pip install -e ./tropic-core -e ./tropic-api[strict] --system

EXPOSE 8000
CMD ["tropic-api"]
