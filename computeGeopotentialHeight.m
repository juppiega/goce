function [addStruct] = computeGeopotentialHeight(addStruct)
% TODO: Lisaa leveyspiirin ja pyorimisen vaikutus.

R = 6356770;
z0 = 130E3;
addStruct.Z = (R + z0) * (addStruct.altitude - z0 * 1E-3) ./ (R + addStruct.altitude*1000);

end