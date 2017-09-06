omit = --omit '*indigo*'
coverage = coverage run $(omit) --source antismash -m pytest
integration_flags = --override-ini=python_files=integration_*.py
integration_coverage = .coverage_integration
sanity_run = echo "sanity TTA run" && ./run_antismash.py --minimal ../inputs/nisin.gbk --tta

unit:
	echo "sanity TTA run" && ./run_antismash.py --tta antismash/test/integration/data/nisin.gbk
	pytest antismash

integration:
	python -m pytest antismash $(integration_flags)

clean:
	find . -name '*.pyc' | xargs rm
	find . -name '__pycache__' | xargs rm -rf

cover: coverage

combined-coverage: coverage
	COVERAGE_FILE=$(integration_coverage) $(coverage) $(integration_flags)
	coverage combine .coverage $(integration_coverage)
	coverage html -d cover
	coverage report

coverage:
	$(sanity_run)
	rm -rf cover .coverage $(integration_coverage)
	$(coverage) antismash
	coverage html -d cover
	coverage report 

