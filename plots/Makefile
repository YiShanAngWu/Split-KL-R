SHELL=/bin/bash
.PHONY: rho

plot:
	python plot.py
	#cd table; pdflatex All_bound_table_pdf.tex
	#cd table; pdflatex All_Lntest_table_pdf.tex
	#cd table; pdflatex All_bestSigma2_table_pdf.tex
	#cd table; pdflatex FW_FWEL_table_pdf.tex
	#cd table; pdflatex Vanilla_FW_table_pdf.tex
	#cd table; pdflatex FWEL_values_pdf.tex
	cd figure; pdflatex multi_graph.tex
	cd table; rm -f *.{log,aux}
	cd figure; rm -f *.{log,aux}
	rm -f texput.log

clean:
	cd table; rm -f *.pdf
