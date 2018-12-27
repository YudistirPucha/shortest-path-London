close all ;clear all;clc;
[t, stations]= xlsread('Final Data.xlsx',1,'A2:A303');%Station Names 
stationplatforms=xlsread('Final Data.xlsx',1,'B2:B303'); % Number of platforms at each station 
network=xlsread('Final Data.xlsx',2,'A2:D407'); % src, dest, line, travel time
linetransfers = xlsread('Final Data.xlsx',3, 'A2: C208'); %Interchanging times 
gatetoplatform = xlsread('Final Data.xlsx',5, 'A2: C416'); 
[~,tubeLines] = xlsread('Final Data.xlsx',4,'B2:B14');


% Go through every edge in the network and add 2 corresponding platforms
% Leave only unique ones and compare against the actual number of platforms
platforms=zeros(1,2*size(network,1));
for i=1:size(network,1)
    platforms(2*i-1) = network(i,1) + 1000*network(i,3);
    platforms(2*i) = network(i,2) + 1000*network(i,3);
end
platforms=unique(platforms);

% Now that we know how many platforms are there, we can start building the adjacency matrix
N = length(stations)+length(platforms); % the total number of nodes in the network
nodes=[1:length(stations) platforms]; % these are gates + platforms
A = zeros(N,N); % our adjacency matrix
% First, we need to add actual tube travel edges
for i=1:size(network,1)
    p1=network(i,1)+1000*network(i,3); i1=find(nodes==p1);
    p2=network(i,2)+1000*network(i,3); i2=find(nodes==p2);
    A(i1,i2)=network(i,4); A(i2,i1)=network(i,4);
end

%Third, we need to add transfers between platforms and gates
%adding travel times from gate to each platfrom at station
for j =1:length(gatetoplatform(:,1))
    j2 = find(nodes == gatetoplatform(j,2));
    A(gatetoplatform(j,1),j2) = gatetoplatform(j, 3);
    A(j2 , gatetoplatform(j,1)) = gatetoplatform(j, 3);
end 

% Add line transfers
for  k = 1 :length(linetransfers)
     k1 = find(nodes == linetransfers(k,1));
     k2 = find(nodes == linetransfers(k,2));
     
     A(k1,k2) = linetransfers(k,3);
     A(k2,k1) = linetransfers(k,3);

 end 
A (A==0) = NaN;%Replace 0 with NaN

% Dialog boxes to ask about starting station and end station
d = dir;
str = stations;
[p,~] = listdlg('PromptString','Select starting station:',...
                'SelectionMode','single',...
                'ListString',str);
            
k = dir;
str = stations;
[q,~] = listdlg('PromptString','Select end station:',...
                'SelectionMode','single',...
                'ListString',str);
%Initialising Dijsktra's Algorithm          
%Set 

for i = 1:max(size(nodes))
    dist(i) = inf;
    % set c is 3 set b is 2 and set a is 1
    label(i) = 3;
    prev(i) = 0;
end
    dist(p) = 0;
    label(p) = 2 ;    
%  Dijktra's Algorithm 
while true
  %Finding the mininum distance 
  setB =  min(dist(label ==2));
  %If set B is empty then break out of while loop
  if length(setB)< 1
      break;
  end
      
   znode= find(dist == setB);
   fprintf('%i, znode');
   
   %Randomly select a node with equal distance from current vertex to
   %adjacent vertex
   if length(znode) > 1
       znode = znode(randi(numel(znode)));
   end
   %Finding the adjacent connected to current node znodes
   ynodes = find((A(znode,:) >0) & label > 1);
  
   %Calculating the shortest path for each adjacent vertex 
   for j = 1:length(ynodes)
    
        if label(ynodes(j)) == 3
             dist(ynodes(j)) = dist(znode) + A(znode, ynodes(j));
                label(ynodes(j))= 2;
                  prev(ynodes(j))= znode;
        else
            %check for better path
           	 if (label(ynodes(j))== 2 && dist(ynodes(j)) > dist(znode)+A(znode,ynodes(j)))
                dist(ynodes(j)) = dist(znode) + A(znode,ynodes(j));
                prev(ynodes(j)) = znode;
             end
        end 
   end 
 % Move current vertex to set A 
   label(znode) = 1;
 %Stop when end station is met
   if znode == q
        break;
   end 
  
end

%Display travel time in  dialog box to show message  
s = strcat('The fastest time from  ', {' '} , stations(p), {' '} ,'to  ' , {' '}, stations(znode), ' is ' ,{' '}, num2str(dist(q)), ' minutes' );
h=msgbox (s , 'The Time');

%Finding the path  
w =[];
k = 1;
w(k) = q;

%Collecting the routes   
 for k = 2: max(prev)
  
    q = prev(q);
    w(k) = q ;
    
    if q == 0 
        break 
    end 
   
 end
 
 %convert path from int to string 
 w = fliplr(w);
 v = w >0;
 v = w(v);
 % How to read the last three digits of the sations to identify the station
 % names 
 d =1000;
 v =nodes(v);
 v1 = mod(v,d);
 linesNumbers = (v - v1) / d;
 linesNumbers = unique(linesNumbers,'stable');

% Take out the zeros 
 linesNumbers(linesNumbers ==0) = []; 
 Lines = tubeLines(linesNumbers);
 Lines = insertBefore(Lines, 1 , ' take');  % works only with 2017 version
%  of MATLAB

 gatetoplatform = (v1 - v) / d;
 gatetoplatform = unique(gatetoplatform);
 
 %Show interchanging stations 
 v1 =  v1(diff(v1)==0); 
 u = stations(v1);

% InsertBefore "Take" before lines
 u1 = strjoin(u , Lines);
 u1 = insertAfter(u1, 'Line' , ' to '); % works only with 2017 version
%  of MATLAB

h=msgbox((u1 ), 'The route');
