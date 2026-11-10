SHELL := /bin/bash

help:
	@echo "Makefile for epiworld models"
	@echo ""
	@echo "Available targets:"
	@echo "  all         - Render all diagrams"

all:
	@for file in *.mmd; do \
		echo "Processing $$file..."; \
		curl -sSf --retry 5 --retry-delay 3 -o "$${file%.mmd}.png" \
			"https://mermaid.ink/img/$$(base64 < "$$file" | tr -d '\n')" || \
			{ rm -f "$${file%.mmd}.png"; echo "Failed to render $$file"; exit 1; }; \
		sleep 2; \
	done && \
	rsync -avz *png ../../../docs/assets/img/ && \
	rm -f *.png
