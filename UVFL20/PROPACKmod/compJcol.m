%%***************************************************
%% find the number of nonzero elements (cumulatively) 
%% in each column of a sparse matrix
%%
%% NNLS, version 0: 
%% Copyright (c) 2009 by
%% Kim-Chuan Toh and Sangwoon Yun 
%%***************************************************

   function jcB = compJcol(JJ)

   jcBtmp = [0; find(diff(JJ))];
   
   coltmp = JJ(jcBtmp+1); 
   
   cnt = 0;  
     
   for k=1:length(jcBtmp) 
       
       if (k==1) 
          len = coltmp(k);
       else
          len = coltmp(k)-coltmp(k-1); 
       end
       
       jcB(cnt+[1:len]) = jcBtmp(k)*ones(1,len);
      
       cnt = cnt + len;      
   end
          
%      mm1 = floor(length(jcBtmp)/4);
%      
%      for k=1:2
%          
%          len = coltmp(k); 
%              
%          jcB(cnt+[1:len]) = jcBtmp(k)*ones(1,len);
% 
%          cnt = cnt + len;
%      end
%      
%      for k = 3:mm1
%          
%            if k==3
%               cnt = cnt-1;
%            end
%            
%            len = coltmp(k)-coltmp(k-2);
%              
%            jcB(cnt+[1:len]) = jcBtmp(k)*ones(1,len);
%       
%            cnt = cnt + len;
%      end
%      
%      for k = mm1+1:length(jcBtmp)
%          
%          if k==mm1+1
%          
%              cnt = cnt-2;
%          
%          end
%          
%          len = coltmp(k);
%              
%          jcB(cnt+[1:len]) = jcBtmp(k)*ones(1,len);
%       
%          cnt = cnt + len;
%          if cnt>length(JJ)
%              break;
%          end    
%      end      
       
    jcB = [jcB,length(JJ)];
 
%%***************************************************
