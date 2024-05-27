function output = circmask_part(r,x1,x2,s1,s2)

% Binary mask with 1 for any pixel with distance from center < r+sqrt(0.5). 
% Roughly matches ImageJ circmask definition.
% Circmask is created around coordinates [x1,x2] in a full frame of size [s1,s2].
% __________________________________
%
% copyright: Fiona Muellner, Institute of Molecular and Clinical Ophthalmology Basel, 24.5.2024


output = zeros(s1,s2,'logical');
for i = 1:s1
    for j = 1:s2
        if sqrt((x1-i)^2 + (x2-j)^2) < r+sqrt(0.5)
            output(i,j) = 1;
        end
    end
end
