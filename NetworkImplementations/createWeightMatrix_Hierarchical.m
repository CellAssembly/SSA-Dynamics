% script to create functional connectivity matrices hierarchical

Eneuron = EneuronNum;
Ineuron = IneuronNum;

num_groups_l1 = 16;
subgroups = 2;
num_groups_total = subgroups*num_groups_l1; 

Eneuron_l1 = Eneuron/num_groups_l1; % make sure this divides properly
Eneuron_l2 = Eneuron_l1/subgroups;


%%%%%%%%
% Probabilities

pEI =.5;
pIE =.5;
pII =.5;

% EE couplings
% pEE= p1/REE;      %.182
% pEE_1 = p1/REE;     %.2
% pEE_2 = p1;    %.22


REE_l1 = 1.45;%1.7;
REE_l2 = 3.7;%3.2;

p1 = .2/ ( 1/num_groups_total + (subgroups-1)/(num_groups_total*REE_l2) + (1-1/num_groups_l1)/(REE_l1*REE_l2) );


pEE= p1/(REE_l1*REE_l2);  
pEE_1 = p1/REE_l2;  
pEE_2 = p1; 

%%%%%%%%%
% Weights
% 

w_EI     = 1/sqrt(mult)*0.042; 
w_IE     = 1/sqrt(mult)*0.0105;
w_II     = 1/sqrt(mult)*.042;

%% EE weights
wEE = 1/sqrt(mult)*.022;   
wEE_1 = 1/sqrt(mult)*.022;
wEE_2 = 1/sqrt(mult)*.023; 


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% procedure for II IE and EI w/o hierarchies
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% AIE = (rand(Ineuron,Eneuron_l1)<pIE);
% AIE = repmat(AIE,1,num_groups_l1);
AIE = (rand(Ineuron,Eneuron)<pIE);
weightsIE = AIE.*(ones(size(AIE)).*w_IE);

% AEI = (rand(Eneuron_l1,Ineuron)<pEI);
% AEI = repmat(AEI,num_groups_l1,1);
AEI = (rand(Eneuron,Ineuron)<pEI);
weightsEI = AEI.*(ones(size(AEI)).*w_EI);


AII = (rand(Ineuron)<pII);
weightsII = AII.*w_II;


%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% CASE PERFECT Symmetry
%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%
% Part 1 -- connectivity


s= cell(1,num_groups_l1);
% comment out if the same grouping everywhere...
for i=1: num_groups_l1
    
    
    ss = cell(1,subgroups);
    A_sub = zeros(Eneuron_l1);
    A_sub2 = zeros(Eneuron_l1);
    
    for ii= 1:subgroups
        % subgroup
        A_diag = (rand(Eneuron_l2)<pEE_2);
        A_diag(A_diag>0) = 1;         % make sure we just encode topology, no effect if unsymmetric version
        ss{ii}= A_diag;
    end
    
    % comment out to get exactly the same
    % assemble cross couplings of matrix
    for ii=1:subgroups-1
        for j=ii:subgroups-1
            
            % cross groups
            A_offdiag1 =(rand(Eneuron_l2)<pEE_1);
            A_offdiag1(A_offdiag1>0) = 1;
            % A_offdiag1 = (A_offdiag1 + A_offdiag1')/2; % symmetrise optional?!
            A_offdiag1(A_offdiag1>0) = 1;         % make sure we just encode topology
            
            % comment out to get exactly the same
            % assemble cross couplings of matrix
            % for ii=1:num_groups_l1-1
            %     for j=i:num_groups_l1-1
            A_sub(Eneuron_l2*(ii-1)+1:Eneuron_l2*ii ,Eneuron_l2*j+1:Eneuron_l2*(j+1)) = ...
                A_offdiag1;
        end
    end
    
        % comment out to get exactly the same
    % assemble cross couplings of matrix
    for ii=1:subgroups-1
        for j=ii:subgroups-1
            
            % cross groups
            A_offdiag1 =(rand(Eneuron_l2)<pEE_1);
            A_offdiag1(A_offdiag1>0) = 1;
            % A_offdiag1 = (A_offdiag1 + A_offdiag1')/2; % symmetrise optional?!
            A_offdiag1(A_offdiag1>0) = 1;         % make sure we just encode topology
            
            % comment out to get exactly the same
            % assemble cross couplings of matrix
            % for ii=1:num_groups_l1-1
            %     for j=i:num_groups_l1-1
            A_sub2(Eneuron_l2*(ii-1)+1:Eneuron_l2*ii ,Eneuron_l2*j+1:Eneuron_l2*(j+1)) = ...
                A_offdiag1;
        end
    end
    
    
    A_sub = A_sub+ A_sub2' + blkdiag(ss{:});
    
    
    % comment out if the same grouping everywhere...
    % s= cell(1,num_groups_l1);
    % for i=1: num_groups_l1
    
    s{i} = A_sub;
end

A = zeros(Eneuron);
A2 = A;
for i=1:num_groups_l1-1
    for j=i:num_groups_l1-1  

% cross groups
A_cross = (rand(Eneuron_l1)<pEE);
% A_cross = (A_cross + A_cross')/2;
A_cross(A_cross>0) = 1;

% comment out to get exactly the same
% assemble cross couplings of matrix
% for i=1:num_groups_l1-1
%     for j=i:num_groups_l1-1       
        A(Eneuron_l1*(i-1)+1:Eneuron_l1*i ,Eneuron_l1*j+1:Eneuron_l1*(j+1)) = ...
            A_cross;    
        
        
% cross groups2
A_cross = (rand(Eneuron_l1)<pEE);
% A_cross = (A_cross + A_cross')/2;
A_cross(A_cross>0) = 1;

% comment out to get exactly the same
% assemble cross couplings of matrix
% for i=1:num_groups_l1-1
%     for j=i:num_groups_l1-1       
        A2(Eneuron_l1*(i-1)+1:Eneuron_l1*i ,Eneuron_l1*j+1:Eneuron_l1*(j+1)) = ...
            A_cross;  
    end
end

% old symmetric version
% A = A+A'+ blkdiag(s{:});

A = A+A2'+ blkdiag(s{:});

%%%%%%%%%%%%%%%%%%
% Part 2 -- just weights

% two level indicator matrices
pvectorl1 = kron(1:num_groups_l1,ones(1,Eneuron_l1));
pvectorl2 = kron(1:num_groups_total,ones(1,Eneuron_l2));

H1 = transformPartitionVectorToHMatrix(pvectorl1);
H2 = transformPartitionVectorToHMatrix(pvectorl2);

W_2 = (ones(Eneuron)*wEE_2 ) .* (H2*H2');
% example random version below...
%W_2 = (ones(Eneuron)*wEE_2*(0.5 + rand(Eneuron))) .* (H2*H2');

W_1 = (ones(Eneuron)*(wEE_1) ) .* (H1*H1');
W_1(logical(H2*H2')) = 0;


W = ones(Eneuron)*(wEE);
W(logical(H1*H1')) = 0;

% total weight matrix
W = W+W_1+W_2;

%%%%%%%%%%%%%%%%%%%
% Final part combine weights with topology

weightsEE = A.*W;
