function T_N = divMat(obj)
%DIVMAT Solves system of equations using LU factors.
T_N =  obj.luFac.Q*(obj.luFac.U\(obj.luFac.L\(obj.luFac.P*(obj.luFac.R\obj.D))));
end

