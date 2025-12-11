SHELL := /bin/bash

help:
	@echo "Makefile for epiworld models"
	@echo ""
	@echo "Available targets:"
	@echo "  all         - Render all diagrams"

all:
	@for file in *.mmd; do \
		echo "Processing $$file..."; \
		cat "$$file" | base64 -w0 | xargs -I {} curl -s "https://mermaid.ink/img/{}" > "$${file%.mmd}.png"; \
	done
