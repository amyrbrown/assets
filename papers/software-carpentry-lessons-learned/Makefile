CONVERT=pandoc
SCLL=software-carpentry-lessons-learned

all : ${SCLL}.html

${SCLL}.html : ${SCLL}.md
	${CONVERT} < $< > $@

pdf : ${SCLL}.md
	@echo '\\documentclass{article}' > ${SCLL}.tex
	@echo '\\usepackage{graphicx}' >> ${SCLL}.tex
	@echo '\\usepackage{hyperref}' >> ${SCLL}.tex
	@echo '\\begin{document}' >> ${SCLL}.tex
	${CONVERT} -t latex $< | sed -e 's/begin{enumerate}\[1.\]/begin{enumerate}/' >> ${SCLL}.tex
	@echo '\\end{document}' >> ${SCLL}.tex
	pdflatex ${SCLL}.tex
	pdflatex ${SCLL}.tex

clean :
	@rm -f ${SCLL}.html ${SCLL}.aux ${SCLL}.log ${SCLL}.out ${SCLL}.pdf ${SCLL}.tex *~
