function logRes = structureComp(a,b,inds)

% Check if two structures are exactly the same (i.e. same fieldnames and
% values within each fieldname
%
% Usage: logRes = structureComp(a,b)

%% First check that each struct has the same number of fieldnames

if isempty(inds)
    
    logRes = isequaln(a,b);

end

%% Next check if each fieldname holds the same value

if ~isempty(inds)
   
    % Get a list of all the structure fieldnames
    fnList = fieldnames(a);
    
    % If you want to just check a few fields, provide their indices 
    if ~isempty(inds)
       fnList = fnList(inds); 
    end
    
    for i = 1:numel(fnList)
       
        % Special handling if any of the fields are cell arrays
        % (will fail if cell array is full of cell arrays... don't do that)
        if iscell(a.(fnList{i}))
            
            atemp = a.(fnList{i});
            btemp = b.(fnList{i});
            
            for j = 1:numel(atemp)
               
                temp(j) = all(all(atemp{j} == btemp{j}));
                
            end
            
            test(i) = all(temp);
            
        elseif isstruct(a.(fnList{i}))
            
            atemp = a.(fnList{i});
            btemp = b.(fnList{i});
            
            %%%% is this good practice to call recursively?
            test(i) = structureComp(atemp,btemp,[]);
            
        else
            % Can just do other arrays as a big logical block with any
            test(i) = all(a.(fnList{i}) == b.(fnList{i}));
        end
        
    end
   
    % If any of the fieldnames are different, note it
    if all(test)
        logRes = true;
    else
        logRes = false;
    end
    
end

end