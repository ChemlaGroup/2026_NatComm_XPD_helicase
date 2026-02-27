function [stgTime, stgPos] = StageMoveTime(Date, stgTimeFilename)

% Read _stagetime.txt files generated during time and initialize start of
% save time to 0 and outputs an array. 
% _stagetime.txt saves start of data save, starts/ends moves and end of
% data save. 
% Created 02/10/19 SY 

global datadirectories


% FileStringArr = strsplit(stgTimeFilename, '_');
% Date = FileStringArr(1);
% stgTimeFilePath = [datadirectory Date '/' stgTimeFilename]

stgTimeFilePath = [datadirectories Date '/' stgTimeFilename];

timeArr = textread(stgTimeFilePath, '%s');

stgTime =[];
stgPos = ['start'];
for i = 1:length(timeArr) 

%        disp(timeArr{i})
   timestrarr = strsplit(timeArr{i},':'); 
   if length(timestrarr) > 1 
       timeHr = 60*60*str2num(timestrarr{1}); 
       timeMin = 60*str2num(timestrarr{2}); 
       timeSec = str2num(timestrarr{3}); 
       currtime = timeHr + timeMin + timeSec;
       stgTime = [stgTime; currtime]; 
   else
       stgPos = [stgPos ' ' timestrarr{1}] ;
   end
end
stgTime = stgTime - stgTime(1); 
%stgPos = [stgPos ' ' 'end'];
stgPos = strsplit(stgPos, ' '); 
stgPos = stgPos(2:end);

end