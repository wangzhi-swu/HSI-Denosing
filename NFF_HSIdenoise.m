% This program is free software: you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation, either version 3 of the License, or
% (at your option) any later version.
%
% This program is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
%
% You should have received a copy of the GNU General Public License
% along with this program.  If not, see <http://www.gnu.org/licenses/>.
%
% -------------------------------
%   version 1.0, August 1, 2016
% -------------------------------

% This program was based on the code
% Copyright(c) <Yongyong Chen> version 1.0,  March, 23, 2017. 
% -------------------------------
%   version 1.0, Oct 26, 2022
% -------------------------------

function  [ output_image ] = NFF_HSIdenoise( oriData3_noise,M,stepszie,a,setmu,setlambda)

% solver for Nonconvex fraction function
%
%        min  P(a) + \lambda ||S||_1, 
%        st.    ||D-L-S||<= \eta
[m,n,p] = size(oriData3_noise);
clear_image=zeros(m,n,p);
opts.type = 1;
opts.rate = 1.5; % [default 1.1, set [1.1,1.5] to speedup]
opts.tol = 1e-7; % [set 1e-5 to speedup] 
opts.max_itr = 1000;
opts.sig = zeros(min(M^2,p),1);
opts.mu = setmu;
opts.lambda = setlambda;

   
Weight = zeros(m,n,p);      % weight matrix 
patch_block = zeros(M^2,p);  % the size of each block    
R         =   m-M+1; % final row start
C         =   n-M+1; % final column start
rr        =   [1:stepszie:R];% each patch start row position
rr        =   [rr rr(end)+1:R];
cc        =   [1:stepszie:C];
cc        =   [cc cc(end)+1:C];% each patch start column position
row       =   length(rr);
column    =   length(cc);


for   rownumber =1:row
    for columnnumber = 1:column
         i = rr(rownumber);
         j = cc(columnnumber); 
        for  k=1:1:p                     
         patch_reference = oriData3_noise(i:i+M-1,j:j+M-1,k); 
         patch_block(:,k) =  patch_reference(:);
         Weight(i:1:i+M-1,j:1:j+M-1,k) = Weight(i:1:i+M-1,j:1:j+M-1,k)+1;               
        end
        tic;
        [clear_patch_block,S] = NFF_ALM(patch_block,opts,a);
        toc;
        for m2=1:1:p
          clear_image(i:1:i+M-1,j:1:j+M-1,m2) = reshape(clear_patch_block(:,m2),M,M)+clear_image(i:1:i+M-1,j:1:j+M-1,m2);
        end 
     end            
end  
 Weight_last = 1./Weight;
 output_image = Weight_last.*clear_image;
end
