malaria_serial_interval
=======================

This repository contains code for generating Figures 1-5 and implementing an interactive calculator for generation and serial intervals from the following manuscript:

Huber JH, Johnston GL, Greenhouse B, Smith DL, Perkins TA. (2016) **Quantitative, model-based estimates of variability in the generation and serial intervals of _Plasmodium falciparum_ malaria**. *Malaria Journal* 15:490. doi:[http://malariajournal.biomedcentral.com/articles/10.1186/s12936-016-1537-6](10.1186/s12936-016-1537-)

All code contained within this repository is released under the [CRAPL v0.1 License](http://matt.might.net/articles/crapl/).

================

### Figure 1

Within the folder ./code/R/, first execute `source('run_Fig1.R')` and then `source('fig_Fig1.R')` in **R** to produce ./output/fig1.tiff.

### Figure 2

Within the folder ./code/R/, first execute `source('run_Fig2.R')` and then `source('fig_Fig2.R')` in **R** to produce ./output/fig2.tiff.

### Figure 3

Within the folder ./code/R/, first execute `source('run_Fig3.R')` and then `source('fig_Fig3.R')` in **R** to produce ./output/fig3.pdf.

### Figure 4

Within the folder ./code/R/, first execute `source('run_Fig4.R')` and then `source('fig_Fig4.R')` in **R** to produce ./output/fig4.pdf.

### Figure 5

Within the folder ./code/R/, first execute `source('run_Fig5.R')` and then `source('fig_Fig5.R')` in **R** to produce ./output/fig5.pdf.

### SI.calculator

Within the folder ./code/R/, execute `source('SI.calculator.R')` in **R**. Calls to the function `SI.calculator(casehistory, output.type, SI.boolean)` should be made from the command line of **R**. The parameter 'casehistory' is a vector of integer values 1 and 2, where 1 represents a treated or symptomatic case and 2 represents an untreated or asymptomatic case. The parameter 'output.type' is a string that can be set to "PDF", "Gamma" or "NB". If set to "PDF", `SI.calculator()` returns a data frame composed of 't' and 'p', which represent timing and probability densities. If set to "Gamma", `SI.calculator()` returns the fitted shape and scale parameters for the gamma distribution. If set to "NB", `SI.calculator()` returns the fitted size and probability parameters for the negative binomial distribution. The parameter 'SI.boolean' takes boolean values. If set to 'TRUE', the function returns the distribution or fitted parameters for the serial interval. If set to 'FALSE', the function returns the distribution or fitted parameters for the generation time interval. 
