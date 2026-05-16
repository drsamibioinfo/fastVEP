# Build stage
FROM rust:latest AS builder

# Install build dependencies
RUN apt-get update && apt-get install -y \
    pkg-config \
    libssl-dev \
    zlib1g-dev \
    && rm -rf /var/lib/apt/lists/*

# Create a new empty shell project
WORKDIR /usr/src/fastvep

# Copy the entire workspace
COPY . .

# Build only CLI binary in release mode
RUN cargo build --release -p fastvep-cli

# Runtime stage
FROM debian:bookworm-slim

# Install runtime dependencies
RUN apt-get update && apt-get install -y \
    libssl3 \
    zlib1g \
    ca-certificates \
    && rm -rf /var/lib/apt/lists/*

# Create application directories
RUN mkdir -p /data /genomes /sa_databases

# Copy binary from the builder stage
COPY --from=builder /usr/src/fastvep/target/release/fastvep /usr/local/bin/

# Set environment variables with defaults
ENV FASTVEP_DATA_DIR=/genomes
ENV FASTVEP_SA_DIR=/sa_databases

# Default command: run the CLI tool
ENTRYPOINT ["fastvep"]
