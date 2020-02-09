%Ethan Green
%February 9th, 2020
%Remove compartmentalization from model.mets
%in order to process using python for webscraping
%of charge data from BiGG
mets = cell(1,length(model.mets));
for i = 1:length(model.mets)
    str = split(model.mets(i),"[");
    mets(i) = str(1);
end
