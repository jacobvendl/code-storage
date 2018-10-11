%this function reads input from a properly formatted .txt file
function [numFM, coordF, magF, coordM, magM, coordS, typeF, typeM, numofMoments] = readInput2(input)

% Open the input file
fileID = fopen(input, 'r');

% read in textline and ignore
textline = fgets(fileID);

% read in number of external forces and moments
textline = fgets(fileID);
ReadFM = str2num(textline);
numFM(1,:) = ReadFM;

% read in textline and ignore
textline = fgets(fileID);
% read in textline and ignore
textline = fgets(fileID);

% read in coordinates of external forces
for i = 1:1:numFM(1,1);
    textline = fgets(fileID);
    ReadPoints = str2num(textline);
    coordF(i,:) = ReadPoints;
end

% read in textline and ignore
textline = fgets(fileID);
% read in textline and ignore
textline = fgets(fileID);

% read in magnitude and direction of external forces
for i = 1:1:numFM(1,1);
    textline = fgets(fileID);
    ReadMags = str2num(textline);
    magF(i,:) = ReadMags;
end

% read in textline and ignore
textline = fgets(fileID);
% read in textline and ignore
textline = fgets(fileID);

% read in location at which external coupled moments are applied
for i = 1:1:numFM(1,2);
    textline = fgets(fileID);
    ReadLocs = str2num(textline);
    coordM(i,:) = ReadLocs;
end

% read in textline and ignore
textline = fgets(fileID);
% read in textline and ignore
textline = fgets(fileID);

% read in magnitude and direction of external coupled moments
for i = 1:1:numFM(1,2);
    textline = fgets(fileID);
    ReadDir = str2num(textline);
    magM(i,:) = ReadDir;
end

% read in textline and ignore
textline = fgets(fileID);
% read in textline and ignore
textline = fgets(fileID);

% read in location of supports
for i = 1:1:6;
    textline = fgets(fileID);
    ReadSup = str2num(textline);
    coordS(i,:) = ReadSup;
end

% read in textline and ignore
textline = fgets(fileID);
% read in textline and ignore
textline = fgets(fileID);



% read in type and direction of reaction

counterF = 1;
counterM = 1;
typeM=[0 0 1];
for i = 1:1:6;
    textline = fgets(fileID);
    if textline(1) == 'F'
        textline(1)=[];
        ReadType = str2num(textline);
        typeF(counterF,:) = ReadType;
        counterF = counterF + 1;
    end
    if textline(1) == 'M'
        textline(1)=[];
        ReadType = str2num(textline);
        typeM(counterM,:) = ReadType;
        counterM = counterM + 1;
    end  
end
numofMoments = counterM -1;