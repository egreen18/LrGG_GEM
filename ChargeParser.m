%Ethan Green
%February 9th, 2020
%Process raw charge data scraped from the web 
%to be a vector of numbers suitable for insertion
%into model
%% Initial parsing of html formatted data
clc
clear
load ChargeRaw.mat
dim = size(charge);
chargestr = cell(1,dim(1));
for i = 1:dim(1)
    str = split(charge(i,:),[">","<"]);
    chargestr{i} = str{3};
end
clear i str
%% Identifying metabolites with mutliple charges and asking for user input
load MetParsed.mat
load LrGG_Model.mat
if isfile('ChargeInt.mat') == 1
    load ChargeInt.mat
    disp("Loaded previous data, please continue from where you left off.")
else
    disp("No previous save data located, initizalizing new save.")
end
URLbase = "http://bigg.ucsd.edu/universal/metabolites/";
start = 1;
for i = start:length(chargestr)
    if length(chargestr{i})>2
        URL = strcat(URLbase,mets(i));
        web(URL);
        str = split(chargestr{i},',');
        disp("PROGRESS: "+num2str(i)+"/"+num2str(length(chargestr)))
        disp("Metabolite "+mets(i)+" has multiple charges: "+str(1)+" and "+str(2)+".")
        disp("It exists in the model as "+model.metFormulas(i)+".")
        newcharge = input("Please choose the appropriate charge:");
        chargestr{i} = num2str(newcharge);
        start = i;
        save ChargeInt.mat %Saving progress. Please push this to git if you make valuable progress.
        clc
    end
end
clear i newcharge
%% Exporting charge data to model
perm = input("Would you like to export the new charge data? [1,0]");
if perm == 1
    model.metCharges = zeros(1,length(chargestr));
    for i = 1:length(chargestr)
        model.metCharges = str2double(chargestr{i});
    end
    disp("Exported. Check for consistency and save file.")
elseif perm == 0
    disp("Rerun section when ready to export data.")
else
    disp("Invalid response type. Please rerun section.")
end
clear i perm