%
%
% File author(s):  <Julian Andres Garcia Grajales>, (C) 2014
%
% Copyright: this software is licenced under the terms stipulated in the license.txt file located in its root directory
%
%


clear all;
close all;

% This is the input for this file
tree = load('neuron_cut_from_neuromorpho.swc');
% Output of this file and input for Neurite
neuron_for_neurite = fopen('neuron_for_neurite.txt','w');

cont=0;
cont2=0;
% This part is to load the file
for i=1:size(tree,1)
   finding=tree(i,7);
   vake=find(tree(:,7)==finding);
   if(size(vake,1)>1)
     cont=cont+1;
     indexing(cont) = vake(2);
     voluta=find(indexing(cont)==indexing(:));
     values(cont)=tree(indexing(cont),7);
     if(size(voluta,1)>1)
         
     else
         cont2=cont2+1;
         indexing_true(cont2)=indexing(cont);
         values_true(cont2)=tree(indexing(cont),7);
     end
     
   end
end

for i=1:size(values_true,2)
real_index(i)=find(tree(:,1)==values_true(i));
end

% Refining the mesh
% Firstly creating the whole tree
for i=1:size(tree,1)

index(i,1) = i; % index
index(i,2) = i-1; % mother
index(i,3) = i+1; % daughter
index(i,4) = 0; % Flag for branching. This value will be overwritten just in case
index(i,5) = 0;  % If I am branching this is my left child.  This value will be overwritten just in case
index(i,6) = tree(i,3); % x
index(i,7) = tree(i,4); % y
index(i,8) = tree(i,5); % z
index(i,9) = tree(i,6); % diameter!!!!!
end

index(1,2)=-1; % To keep the swc format

% Here we modify the values for the special elements:
 for i=1:size(values_true,2)
    index(real_index(i),4)=1; % branching!!!!!!
    index(real_index(i),5)=indexing_true(i); % its left child!!!!
    index(indexing_true(i),2)=real_index(i); % The mother for the left child!!!!
    index(indexing_true(i)-1,3)=indexing_true(i)-1; % for the terminal, my right child is myself. NEURITE condition!!!
 end
 
 
dist_vector_true = zeros(size(tree,1),1);
 for i=2:size(tree,1)
    dist_vector_true(i-1) = distance(index(index(i,2),6),index(index(i,2),7),index(index(i,2),8),index(i,6),index(i,7),index(i,8));
 end

index(1,10)=0;
 
 for i=1:size(tree,1)
index(i,10) = dist_vector_true(i); % diameter!!!!!
end
 % the radius for 1-2 is the 2
 index(1,9) = index(2,9);
 
% Now we define the element size for Neurite
total_elements = size(tree,1)-1 % 

times=0;
numstepsx = total_elements
for i=1:total_elements
    
pos=i; 
    
true_elements(pos,1) = pos; % Me
true_elements(pos,2) = pos-1; % My parent
true_elements(pos,3) = pos+1; % Right child
true_elements(pos,4) = 0; % Branching?
true_elements(pos,5) = 0; % Left child
true_elements(pos,6) = index(i,9); % diameter
true_elements(pos,7) = index(i,10); % lenght;

end


for i=1:total_elements
   
    if(index(i,4)==1)
      pos = i;   
      true_elements(pos,4) = index(i,4); % branching 
      true_elements(pos,5) = (index(i,5)-1)+1; % my left child
      true_elements((index(i,5)-1)+1,2)=index(index(i,5),2);% its parent for left child !!!!        
    end
    if(index(i,1)==index(i,3))
        pos = i; 
       true_elements(pos,3) = true_elements(pos,1);  %terminal elements
        
    end
end
true_elements(numstepsx,3) = true_elements(numstepsx,1);  

if(size(true_elements,1)~=numstepsx)
  error('After the refinement process there is an error. %i must be equal to %i',size(true_elements,1),numstepsx);
end


for i=1:numstepsx
   fprintf(neuron_for_neurite,'%i %i %i %i %i %f %f\n',true_elements(i,1),true_elements(i,2),true_elements(i,3),true_elements(i,4),true_elements(i,5),true_elements(i,6),true_elements(i,7));
    
end


 fclose(neuron_for_neurite);
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 





