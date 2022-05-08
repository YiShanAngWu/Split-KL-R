SHELL=/bin/bash
.PHONY: rho

plot:
	python plot.py
	cd table; pdflatex FW_FWEL_table_pdf.tex
	cd table; pdflatex FWEL_values_pdf.tex
	cd table; rm -f *.{log,aux}
	rm -f texput.log

clean:
	cd table; rm -f *.pdf
	cd table; rm -f *.pdf
