#!/usr/bin/env python3
"""Render the six manuscript figures and refresh their source-data CSVs."""
from __future__ import annotations
import csv, json
from pathlib import Path
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch
ROOT=Path(__file__).resolve().parents[2]; RES=ROOT/'results'; OUT=ROOT/'paper/figures'; SRC=ROOT/'paper/source_data'
OUT.mkdir(parents=True,exist_ok=True); SRC.mkdir(parents=True,exist_ok=True)
def save(fig,name):
    fig.savefig(OUT/f'{name}.png',dpi=200,bbox_inches='tight'); fig.savefig(OUT/f'{name}.svg',bbox_inches='tight'); plt.close(fig)
def rows(path):
    with open(path,encoding='utf-8') as f:return list(csv.DictReader(f))
def write(name,fields,data):
    with open(SRC/name,'w',newline='',encoding='utf-8') as f:w=csv.DictWriter(f,fieldnames=fields);w.writeheader();w.writerows(data)
# Figure 1
wf=[{'order':1,'stage':'Sequence inputs','detail':'Gn, Gc, N\nGenBank-mapped'},{'order':2,'stage':'HLA prediction','detail':'55 class I\n42 class II'},{'order':3,'stage':'Strong-binding filter','detail':'EL percentile\nthresholds'},{'order':4,'stage':'Redundancy collapse','detail':'Overlapping peptides\nto non-redundant cores'},{'order':5,'stage':'Core annotation','detail':'Conservation, HLA breadth,\nhuman-proteome similarity'},{'order':6,'stage':'Fixed candidate set','detail':'15 class I + 10 class II\nrepresentatives'},{'order':7,'stage':'Population analysis','detail':'2,693 genotypes\nbootstrap intervals'}]
write('fig1_workflow.csv',list(wf[0]),wf)
fig,ax=plt.subplots(figsize=(13,3.2));ax.axis('off');w=.125;gap=.012;x=.01
for r in wf:
    box=FancyBboxPatch((x,.25),w,.55,boxstyle='round,pad=0.005,rounding_size=.01',fill=False,linewidth=1.2,transform=ax.transAxes);ax.add_patch(box);ax.text(x+w/2,.61,r['stage'],ha='center',va='center',fontsize=11,fontweight='bold',transform=ax.transAxes);ax.text(x+w/2,.43,r['detail'],ha='center',va='center',fontsize=9,transform=ax.transAxes)
    if r['order']<len(wf): ax.annotate('',xy=(x+w+gap*.65,.525),xytext=(x+w,.525),xycoords=ax.transAxes,arrowprops=dict(arrowstyle='->',lw=1))
    x+=w+gap
save(fig,'fig1_workflow')
# Data
cov=json.loads((RES/'population_coverage.json').read_text()); c1=rows(RES/'final_epitopes_class1.csv'); c2=rows(RES/'final_epitopes_class2.csv'); all1=rows(RES/'cores_class1_scored.csv'); all2=rows(RES/'cores_class2_scored.csv'); h1=json.loads((RES/'human_similarity_class1.json').read_text()); h2=json.loads((RES/'human_similarity_class2.json').read_text())
# Figure 2
stages=['Prediction rows','Unique strong peptides','Non-redundant cores','Selected cores']; valsI=[342430,len([k for k in h1 if k!='_meta']),len(all1),len(c1)]; valsII=[64680,len([k for k in h2 if k!='_meta']),len(all2),len(c2)]
red=[{'class':c,'stage':s,'count':v} for c,vs in [('I',valsI),('II',valsII)] for s,v in zip(stages,vs)]; write('fig2_analysis_reduction.csv',list(red[0]),red)
fig,ax=plt.subplots(figsize=(9,5));x=range(4);ax.plot(x,valsI,marker='o',label='Class I');ax.plot(x,valsII,marker='o',label='Class II');ax.set_yscale('log');ax.set_xticks(list(x));ax.set_xticklabels(stages,rotation=12,ha='right');ax.set_ylabel('Count (log scale)');ax.legend();ax.grid(axis='y',alpha=.25)
for vs in [valsI,valsII]:
    for i,v in enumerate(vs):ax.text(i,v*1.25,f'{v:,}',ha='center',fontsize=9)
save(fig,'fig2_analysis_reduction')
# Figure 3
con=[]
for cls,data in [('I',c1),('II',c2)]:
    for i,r in enumerate(data,1):con.append({'selected_id':f'{cls}{i:02d}','class':cls,'andv_mean_identity_pct':float(r['cons_ANDV_mean_ident'])*100,'full_panel_mean_identity_pct':float(r['cons_all_mean_ident'])*100})
write('fig3_selected_conservation.csv',list(con[0]),con);fig,ax=plt.subplots(figsize=(9,5));labels=[r['selected_id'] for r in con];xs=range(len(con));ax.scatter(xs,[r['andv_mean_identity_pct'] for r in con],label='ANDV-annotated subset');ax.scatter(xs,[r['full_panel_mean_identity_pct'] for r in con],marker='x',label='Full New World panel');ax.axhline(95,ls='--',lw=1);ax.axhline(80,ls=':',lw=1);ax.set_xticks(list(xs));ax.set_xticklabels(labels,rotation=90,fontsize=7);ax.set_ylabel('Mean sequence identity (%)');ax.legend();ax.grid(axis='y',alpha=.2);save(fig,'fig3_selected_conservation')
# Figure 4
br=[]
for cls,data in [('I',c1),('II',c2)]:
    m=cov['representative_allele_breadth']['class_I' if cls=='I' else 'class_II']
    for i,r in enumerate(data,1):br.append({'selected_id':f'{cls}{i:02d}','class':cls,'strong_binding_allele_count':m[r['core_id']]})
write('fig4_representative_allele_breadth.csv',list(br[0]),br);fig,axes=plt.subplots(2,1,figsize=(9,5.2))
for ax,cls in zip(axes,['I','II']):
    d=[r for r in br if r['class']==cls];ax.bar([r['selected_id'] for r in d],[r['strong_binding_allele_count'] for r in d]);ax.set_title(f'Class {cls}',loc='left',fontweight='bold');ax.set_ylabel('Alleles');ax.tick_params(axis='x',rotation=60);ax.grid(axis='y',alpha=.2)
fig.supylabel('Strong-binding HLA alleles per representative');save(fig,'fig4_representative_allele_breadth')
# Figure 5
cv=[]
for s in ['AFR','AMR','EAS','EUR','SAS']:cv.append({'superpopulation':s,'typed_loci_combined_pct':cov['drb1_dqb1']['by_superpopulation_combined'][s]['combined'],'drb1_only_combined_pct':cov['drb1_only']['by_superpopulation_combined'][s]['combined'],'n':cov['drb1_dqb1']['by_superpopulation_combined'][s]['n']})
write('fig5_coverage_by_superpopulation.csv',list(cv[0]),cv);fig,ax=plt.subplots(figsize=(8.5,4.7));x=list(range(len(cv)));w=.36;ax.bar([i-w/2 for i in x],[r['typed_loci_combined_pct'] for r in cv],w,label='DRB1 + DQB1 typed-loci model');ax.bar([i+w/2 for i in x],[r['drb1_only_combined_pct'] for r in cv],w,label='DRB1-only sensitivity model');ax.set_xticks(x);ax.set_xticklabels([r['superpopulation'] for r in cv]);ax.set_ylim(82,100);ax.set_ylabel('Combined genotype coverage (%)');ax.legend();ax.grid(axis='y',alpha=.2);save(fig,'fig5_coverage_by_superpopulation')
# Figure 6
arch=[{'component':'Signal feature','count':1},{'component':'Class-I representatives','count':15},{'component':'Class-II representatives','count':10},{'component':'Predicted B-cell regions','count':2},{'component':'Helper segment','count':1}];write('fig6_candidate_architecture.csv',list(arch[0]),arch);fig,ax=plt.subplots(figsize=(12,2.7));ax.axis('off');labels=['Signal\nfeature','Class-I\nrepresentatives x15','Class-II\nrepresentatives x10','Predicted B-cell\nregions x2','Helper\nsegment'];xs=[.08,.29,.57,.81,.96]
for x,l in zip(xs,labels):ax.text(x,.62,l,ha='center',va='center',fontsize=11,transform=ax.transAxes)
ax.text(.5,.22,'Conceptual organization; complete candidate sequence is distributed separately as FASTA within the research package.',ha='center',fontsize=9,transform=ax.transAxes);save(fig,'fig6_candidate_architecture')
print('wrote 6 figures and 6 source-data CSVs')
