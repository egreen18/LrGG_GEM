%% Loading Data
load genome.mat
load models.mat
%% Pulling protein ID data from genome file
protein_id = zeros(1,length(genome.CDS));
n1 = 0; %Counting the number of proteins without IDs in the genome
for i = 1:length(genome.CDS)
    if isempty(genome.CDS(i).protein_id) %Checking for missing values, ex row 7
        protein_id(i) = "";
        n1 = n1+1;
    else
        data = genome.CDS(i).protein_id;
        str = split(data,["_","."]); %Removing unnecessary information from strings
        protein_id(i) = str2double(str(2));
    end
end
disp(n1+" out of "+i+" proteins are without ID from sequence.gb")
%% Pulling protein ID data from CM model
CMprotein_id = zeros(1,length(modelCM.genes));
n2 = 0; %Counting the number of proteins without IDs in the CM model
for i = 1:length(modelCM.genes)
    str = split(modelCM.genes(i),["_","."]);
    if length(str) < 6 %Checking for completely missing values, ex "spontaneous" on row 28
        CMprotein_id(i) = "";
        n2  = n2+1;
    else
        if contains(str(6),'WP') %Checking for partially missing values, ex row 21
            CMprotein_id(i) = str2double(str(7));
        else
            CMprotein_id(i) = "";
            n2  = n2+1;
        end
    end
end
disp(n2+" out of "+i+" proteins are without ID from the CarveMe model.")
%% Pulling locus tags from genome
locus_id = strings(1, length(genome.CDS));
for i = 1:length(genome.CDS) %Identifying row with /locus_tag
    for j = 1:size(genome.CDS(i).text,1)
        if strfind(genome.CDS(i).text(j,:),'/locus_tag') == 1
            index = j;
            break
        end
    end
    str = split(genome.CDS(i).text(index,:),["/locus_tag=",""""]);
    locus_id(i) = str(3); %Isolating locus tag string
end
disp(i+" locus tags were identified from sequence.gb")
%% Indexing protein ID matches between CM model and genome protein IDs
indexCM = zeros(1,length(CMprotein_id));
n3 = 0; %Counting the number of non-matches between the systems
for i = 1:length(CMprotein_id)
    if isempty(find(protein_id == CMprotein_id(i))>0) %Checking for non-matches
        indexCM(i) = 0;
        n3 = n3+1;
    else
        indexCM(i) = find(protein_id == CMprotein_id(i));
    end
end
disp(n3+" out of "+i+" CM model protein IDs had no match to sequence.gb")
%% Indexing locus ID matches between KB model and genome locus IDs
indexKB = zeros(1,length(modelKB.genes));
n4 = 0; %Counting the number of non-matches between the systems
for i = 1:length(modelKB.genes)
    if isempty(find(locus_id == modelKB.genes(i))>0) %Checking for non-matches
        indexKB(i) = 0;
        n4 = n4+1;
    else
        indexKB(i) = find(locus_id == modelKB.genes(i));
    end
end
disp(n4+" out of "+i+" CM model protein IDs had no match to sequence.gb")
%% Updating CM model with locus tags and protein IDs
new_modelCM = modelCM; 
protein_id = strcat('WP_',strsplit(num2str(protein_id))); %Reconstructing protein ID format
locus_id = cellstr(locus_id); %Normalizing data class
for i = 1:length(indexCM) 
    if indexCM(i) == 0 %Dealing with indexing errors resulting from zeros
        new_modelCM.proteins(i) = {'missing'};
        new_modelCM.genes(i) = {'missing'};
    else
        new_modelCM.proteins(i) = protein_id(indexCM(i)); %Updating values
        new_modelCM.genes(i) = locus_id(indexCM(i));
    end
end
disp("The CarveMe model was updated with gene locus and protein IDs.")
%% Updating KB model with locus tags and protein IDs
new_modelKB = modelKB; 
for i = 1:length(indexKB)
    if indexKB(i) == 0 %Dealing with indexing errors resulting from zeros
        new_modelKB.proteins(i) = {'missing'}; 
    else
        new_modelKB.proteins(i) = protein_id(indexKB(i)); %Updating values
    end
end
disp("The KBase model was updated with protein IDs.")