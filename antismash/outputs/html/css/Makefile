all: bacteria.css fungi.css plants.css

%.css: %.less style.less secmet.less
	lessc --clean-css="--s1 --advanced --compatibility=ie8" $< $@

clean:
	rm *.css

.PHONY: clean
