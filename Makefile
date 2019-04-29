omit = --omit '*/external/*'
coverage = coverage run $(omit) --source antismash -m pytest
integration_flags = --override-ini=python_files=integration_*.py
integration_coverage = .coverage_integration
sanity_run = echo "sanity TTA run" && rm -rf nisin && ./run_antismash.py --minimal antismash/test/integration/data/nisin.gbk

unit:
	$(sanity_run)
	echo "simple reuse test" && ./run_antismash.py --reuse-results nisin/nisin.json
	pytest --durations=3 antismash

integration: clean
	python -m pytest --durations=3 antismash $(integration_flags)

clean:
	rm -f antismash/detection/hmm_detection/data/bgc_seeds.hmm*
	rm -f antismash/modules/clusterblast/data/known/*.dmnd
	rm -f antismash/modules/clusterblast/test/data/*/*.dmnd
	rm -f antismash/outputs/html/css/*.css
	find . -name "*.h3?" -exec rm {} +
	find . -name '*.pyc' | xargs rm -f
	find . -name '__pycache__' | xargs rm -rf

squeakyclean: clean
	find . -name "*.tar.*" -exec rm {} +
	bash -c 'for d in $$(find . -maxdepth 2 -name "index.html"); do DIR=$$(dirname $$d); rm -r $$DIR; done'

cover: coverage

combined-coverage: coverage
	COVERAGE_FILE=$(integration_coverage) $(coverage) $(integration_flags)
	coverage combine .coverage $(integration_coverage)
	coverage html -d cover
	coverage report

coverage:
	$(sanity_run)
	rm -rf cover .coverage $(integration_coverage)
	coverage run $(omit),'*integration_*.py' --source antismash -m pytest antismash
	coverage html -d cover
	coverage report

docker: squeakyclean
	docker build -t antismash/antismash-dev .

.PHONY:	unit integration clean squeakyclean cover coverage combined-coverage docker
