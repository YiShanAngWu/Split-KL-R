SHELL=/bin/bash
.PHONY: rho

plot:
	python plot.py
	cd figure; pdflatex multi_graph.tex
	cd figure; rm -f *.{log,aux}
	rm -f texput.log

clean:
	cd figure; rm -f *.pdf
