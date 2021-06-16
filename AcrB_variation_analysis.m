kpnih=getgenbank('CP008827.1') %get the genome of KPNIH1
dataset=cell(387,80)   %set a dataset
for i=1:387
genome=getgenbank(kpasslist{i})  %Read Klebsiella pneumoniae accession list into the working enviroment first
allCDS=featureparse(genome,'Feature','CDS','Sequence',true)  %featureparse the genome since there will be some errors
seq=cellfun(@nt2aa,{allCDS.Sequence},'uniform',false)
[genome.CDS.translation]=deal(seq{:})
n=numel(genome.CDS)
score=zeros(n,1)
for k=1:n
    score(k)=nwalign(kpnih.CDS(4551).translation,genome.CDS(k).translation)
end   %AcrB is located at 4551 of KPNIH1 genome
dataset{i,1}=kpasslist{i}
a=(score>0)   % finding all AcrB homologes
pos=find(a)   % finding the positions
sc=(score(pos))
new=[pos,sc]
new=sortrows(new,2,'descend')
pos=new(:,1)
sc=new(:,2)
m=length(pos)
for p=1:m
    dataset{i,(p*3)}=sc(p)
    dataset{i,(p*3+1)}=genome.CDS(pos(p)).protein_id
    dataset{i,(p*3+2)}=genome.CDS(pos(p)).translation
end   %write information into the dataset
end
table=cell2table(dataset)
writetable(table,'AcrB_variation.xlsx')