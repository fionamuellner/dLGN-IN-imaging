function output = circmask(r)

% Binary mask with 1 for any pixel with distance from center < r+sqrt(0.5). 
% Roughly matches ImageJ circmask definition.
%
% __________________________________
%
% copyright: Fiona Muellner, Institute of Molecular and Clinical Ophthalmology Basel, 24.5.2024

s = ceil(r);

output = zeros(2*s+1,2*s+1);
for i = 1:2*s+1
    for j = 1:2*s+1
        if sqrt((s+1-i)^2 + (s+1-j)^2) < r+sqrt(0.5)
            output(i,j) = 1;
        end
    end
end
