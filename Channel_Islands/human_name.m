function Out = human_name(S)

% Generates a human-friendlier name for display of species & site names
switch S
    case 'PCLA'
        Out = 'Kelp bass';
    case 'SPUL'
        Out = 'California sheephead';
    case 'SMYS'
        Out = 'Blue rockfish';
    case 'SATR'
        Out = 'Kelp rockfish';
    case 'SMI_TYLER_BIGHT'
        Out = 'Tyler Bight (San Miguel)';
    case 'SMI_CROOK_POINT'
        Out = 'Crook Pt. (San Miguel)';
    case 'SMI_HARRIS_PT_RESERVE'
        Out = 'Harris Pt. Reserve (San Miguel)';
    case 'SMI_CUYLER'
        Out = 'Cuyler (San Miguel)';
    case 'SRI_SOUTH_POINT'
        Out = 'South Point (Santa Rosa)';
    case 'SRI_JOHNSONS_LEE_SOUTH'
        Out = 'Johnsons Lee South (Santa Rosa)';
    case 'SCI_FORNEY'
        Out = 'Forney (Santa Cruz)';
    case 'SCI_PAINTED_CAVE'
        Out = 'Painted Cave (Santa Cruz)';
    case 'SCI_HAZARDS'
        Out = 'Hazards (Santa Cruz)';
    case 'SCI_GULL_ISLE'
        Out = 'Gull Isle (Santa Cruz)';
    case 'SCI_VALLEY'
        Out = 'Valley (Santa Cruz)';
    case 'SCI_YELLOWBANKS'
        Out = 'Yellowbanks (Santa Cruz)';
    case 'SCI_PELICAN'
        Out = 'Pelican (Santa Cruz)';
    case 'SCI_CAVERN_POINT'
        Out = 'Cavern Point (Santa Cruz)';
    case 'SCI_SCORPION'
        Out = 'Scorpion (Santa Cruz)';
    case 'ANACAPA_WEST_ISLE'
        Out = 'West Isle (Anacapa)';
    case 'ANACAPA_MIDDLE_ISLE'
        Out = 'Middle Isle (Anacapa)';
    case 'ANACAPA_EAST_ISLE'
        Out = 'East Isle (Anacapa)';
    case 'ANACAPA_LIGHTHOUSE_REEF'
        Out = 'Lighthouse Reef (Anacapa)';
end