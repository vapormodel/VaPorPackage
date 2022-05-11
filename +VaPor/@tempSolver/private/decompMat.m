function decompMat(obj)
%DECOMPMAT Decomposes T_Solve into LU factors with reordering.
[obj.luFac.L,obj.luFac.U,obj.luFac.P,obj.luFac.Q,obj.luFac.R] = lu(obj.T_Solve);
end

