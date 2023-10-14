function makeHeader(filename, T, R)
%MAKEHEADER Make the extra header file for POV-Ray

if nargin < 3
    R = zeros(1, 3);
end
% Convert to povray coordinates
T(1:2) = -T(1:2);
R(1:2) = -R(1:2);

% Make header
fid = fopen(filename, 'w');
fprintf(fid, '#version version;');
fprintf(fid, '#declare Camera_T = <%.3f, %.3f, %.3f>;\n', T(1)/1000, T(2)/1000, T(3)/1000);
fprintf(fid, '#declare Camera_Rx = %.3f;\n', R(1));
fprintf(fid, '#declare Camera_Ry = %.3f;\n', R(2));
fprintf(fid, '#declare Camera_Rz = %.3f;\n', R(3));
fclose(fid);

end

