import matplotlib
#matplotlib.use('Agg')
import numpy as np
import sys, os
import glob
from shutil import copy, move

def make_web_pages(plotout, htmlout, suites, algo, algo_fname):
    if not os.path.exists(htmlout):
        os.makedirs(htmlout)
    algo_out = os.path.join(htmlout, algo)
    if not os.path.exists(algo_out):
        os.makedirs(algo_out)

    search = os.path.join(plotout, '*'+algo_fname+'*.png')
    print('search ',search)
    plots = glob.glob(os.path.join(plotout, '*'+algo_fname+'*.png'))
    print('len(plots) ',len(plots))
    plots_pdf = glob.glob(os.path.join(plotout, '*'+algo_fname+'*.pdf'))
    print('open ',os.path.join(algo_out, 'assess.html'))
    f = open(os.path.join(algo_out, 'assess.html'), 'w')
    runids = []
    if len(suites) <= 5:
        for run in suites:
            runids.append(run)
    else:
        for run in suites[0:3]:
            runids.append(run)
        runids.append('and others')
    message = 'Assessment of '+' '+' '.join(runids)+' using '+algo_out
    f.write(message+'\n')
    f.write('<head> \n')
    message1 = '<title>'+message+'</title>\n'
    f.write(message1)
    f.write('</head> \n')
    f.write('<body> \n')
    message = '<TABLE cellpadding="4" align="center" border="0" style="border: 1px solid #000000; border-collapse: collapse;">'
    f.write(message+'\n')
    #message = '<center><table border=2 class="mivnchead"><tr> \n'
    #f.write(message)
    # remove old plots first
    for plot in plots:
        print('copy ',plot, os.path.join(htmlout, os.path.basename(plot)))
        move(plot, os.path.join(algo_out, os.path.basename(plot)))
        pname = os.path.basename(plot)
        pfull = os.path.join(algo_out, pname)
        f.write('<TR> \n')
        message = '<TD><a href="'+pname+'"><img src="'+pname+'" alt="Image missing" /></a></TD> \n'
        f.write(message)
        f.write('</TR> \n')
    f.write('</TABLE> /n')
    message= '</td></tr> \n'
    f.write(message)
    f.write('</body> \n')
    f.write('</html> \n')
    f.close()
    for plot in plots_pdf:
        move(plot, os.path.join(algo_out, os.path.basename(plot)))

if __name__ == '__main__':

    suites = ['u-ai674','u-ak681','u-ak687','u-aj059']
    plot_basedir = '/home/users/mjrobert/workspace/storm_analysis_code/plots/'
    figbase = '_'.join([r for r in suites])
    basedir_fig = os.path.join(plot_basedir, figbase)


    htmlout = os.path.join(basedir_fig, 'html')
    if not os.path.exists(htmlout):
        os.makedirs(htmlout)
    
    make_web_pages(basedir_fig, htmlout, suites)
