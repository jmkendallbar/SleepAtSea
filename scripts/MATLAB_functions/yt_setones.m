function res=yt_setones(a,x)
%YT_SETONES outputs start/end indexes of 'one series' in a zeros/ones column vector
%  RES=YT_SETONES(A), for a column vector containing [0 1] values.
%  YT_SETONES (A) returns a two column matrix with the start and end indexes of
%  any 'one serie' in the input vector, respectively.
%  Isolated single ones have the same start and end indexes.
%  YT_SETONES(A,X), returns the same type of index matrix for 'one serie'
%  of size X or more.
%  This function is useful to identify indexes of block of data after
%  using the find command.
%  Example: a=[1;0;0;1;1;0;1;1;1;1;0;0]
%  setones(a)= 1  1
%                         4  5
%                         7  10
%  setones(a,3)= 7  10
%
%  Created by Yann Tremblay on 02 November 2003.
%  Modifications:
% 7 March 2004: treat the case when a is empty
%               replace old NaN outputs by [] outputs.
% 7 March 2004: Bug fixed.
% 24 May 2004: Bug fixed.
% 2 Dec 2004: Bug fixed by M Antolos.


if isempty(a)
    res=[];
    return
else

if nargin==1 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if ~isempty(setdiff(a,[0;1]))
        error('In : setones(a), a must contain zeros or ones ONLY');
    elseif size(a,2)>1
        error('In : setones(a), a must be a column vector');
    end

elseif nargin==2 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if ~isempty(setdiff(a,[0;1]))
        error('In : setones(a,x), a must contain zeros or ones ONLY');
    elseif size(a,2)>1
        error('In : setones(a,x), a must be a column vector');
    elseif ~isnumeric(x)
        error('In : setones(a,x), x must be a numerical value')
    elseif size(x,1)>1 | size(x,2)>1
        error('In : setones(a,x), x must be scalar')
    end
    
elseif nargin>2 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        error('Too many input arguments');
end
    
%%%%% START

    s=size(a,1);
    shift_down=[a(1,1);a(1:s-1,1)];
    difference=a-shift_down;
    borne1=find(difference==1);
    borne2=find(difference==-1)-1;

if isempty(borne1) & sum(a)>1 & a(s,1)==1
    res=[1 , s];
elseif isempty(borne1) & sum(a)>1 & a(s,1)==0
    res=[1 , borne2(1,1)];
elseif isempty(borne1) & sum(a)==0
    res=[];
elseif sum(a)==1
    res=[find(a) , find(a)];  
    
else
    if ~isempty(borne1) & isempty(borne2)
        borne2=s;
    end
    if borne1(1,1)>borne2(1,1)
    borne1=[1;borne1(:,1)];
    end
    if borne1(size(borne1,1),1)>borne2(size(borne2,1),1)
    borne2=[borne2;s];
    end
    res=[borne1 , borne2];
 end
 
if nargin==2 & ~isempty(res)
res=res(find((res(:,2)-res(:,1)+1)>=x),:);  
end
    

end
%%% End of function setones