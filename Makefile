default:
	@echo "--------------------------------------------------------------------"
	@echo "This Makefile supports the following options"
	@echo "--------------------------------------------------------------------"
	@echo "  make install     (Installs JunoPy)"
	@echo "  make install-dev (Installs JunoPy in editable mode)"
	@echo "  make project     (Creates project directory)"
	@echo "  make docs        (Builds the documentation)"
	@echo "  make tests       (Runs all tests in the tests directory)"
	@echo "--------------------------------------------------------------------"

install:
	pip install .

install-dev:
	pip install -e .

project:
	mkdir -p project/data
	mkdir project/spice

docs:
	mkdocs build

tests:
	pytest tests