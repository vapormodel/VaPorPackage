function ret = distPointToLineSegMat(Mat1,Mat2,Pnt)
% Finds shortest distance between a point (pnt) and a line defined
% by two points (vec1 vec2).
l_squared = sum((Mat2 - Mat1).^2,2);

vector1 = (Pnt - Mat1);
vector2 = (Mat2 - Mat1);
dotVector = sum(vector1.*vector2,2);
%T2 = max(0, min(1, dot(Pnt - Vec1, Vec2 - Vec1) / l_squared));
T2 = max(0, min(1, dotVector ./ l_squared));
projection = Mat1 + T2 .* (Mat2 - Mat1);
D2 = vecnorm(Pnt-projection,2,2);
ret = [D2,T2];
end
