OPENQASM 2.0;
include "qelib1.inc";
qreg q[100];
h q[0];
h q[1];
h q[2];
h q[3];
h q[4];
h q[5];
h q[6];
h q[7];
h q[8];
h q[9];
h q[10];
h q[11];
h q[12];
h q[13];
h q[14];
h q[15];
h q[16];
h q[17];
h q[18];
h q[19];
h q[20];
h q[21];
h q[22];
h q[23];
h q[24];
h q[25];
h q[26];
h q[27];
h q[28];
h q[29];
h q[30];
h q[31];
h q[32];
h q[33];
h q[34];
h q[35];
h q[36];
h q[37];
h q[38];
h q[39];
h q[40];
h q[41];
h q[42];
h q[43];
h q[44];
h q[45];
h q[46];
h q[47];
h q[48];
h q[49];
h q[50];
h q[51];
h q[52];
h q[53];
h q[54];
h q[55];
h q[56];
h q[57];
h q[58];
h q[59];
h q[60];
h q[61];
h q[62];
h q[63];
h q[64];
h q[65];
h q[66];
h q[67];
h q[68];
h q[69];
h q[70];
h q[71];
h q[72];
h q[73];
h q[74];
h q[75];
h q[76];
h q[77];
h q[78];
h q[79];
h q[80];
h q[81];
h q[82];
h q[83];
h q[84];
h q[85];
h q[86];
h q[87];
h q[88];
h q[89];
h q[90];
h q[91];
h q[92];
h q[93];
h q[94];
h q[95];
h q[96];
h q[97];
h q[98];
h q[99];
h q[49];
t q[86];
cx q[68], q[87];
h q[58];
cx q[58], q[21];
tdg q[21];
cx q[69], q[21];
t q[21];
cx q[58], q[21];
tdg q[21];
cx q[69], q[21];
t q[21];
cx q[69], q[58];
tdg q[58];
cx q[69], q[58];
t q[69];
t q[58];
h q[58];
h q[67];
cx q[67], q[44];
tdg q[44];
cx q[31], q[44];
t q[44];
cx q[67], q[44];
tdg q[44];
cx q[31], q[44];
t q[44];
cx q[31], q[67];
tdg q[67];
cx q[31], q[67];
t q[31];
t q[67];
h q[67];
h q[98];
t q[32];
h q[16];
cx q[16], q[93];
tdg q[93];
cx q[62], q[93];
t q[93];
cx q[16], q[93];
tdg q[93];
cx q[62], q[93];
t q[93];
cx q[62], q[16];
tdg q[16];
cx q[62], q[16];
t q[62];
t q[16];
h q[16];
h q[35];
h q[60];
h q[88];
cx q[88], q[12];
tdg q[12];
cx q[49], q[12];
t q[12];
cx q[88], q[12];
tdg q[12];
cx q[49], q[12];
t q[12];
cx q[49], q[88];
tdg q[88];
cx q[49], q[88];
t q[49];
t q[88];
h q[88];
s q[43];
h q[45];
cx q[45], q[84];
tdg q[84];
cx q[14], q[84];
t q[84];
cx q[45], q[84];
tdg q[84];
cx q[14], q[84];
t q[84];
cx q[14], q[45];
tdg q[45];
cx q[14], q[45];
t q[14];
t q[45];
h q[45];
h q[85];
h q[39];
cx q[39], q[33];
tdg q[33];
cx q[7], q[33];
t q[33];
cx q[39], q[33];
tdg q[33];
cx q[7], q[33];
t q[33];
cx q[7], q[39];
tdg q[39];
cx q[7], q[39];
t q[7];
t q[39];
h q[39];
h q[55];
cx q[55], q[41];
tdg q[41];
cx q[37], q[41];
t q[41];
cx q[55], q[41];
tdg q[41];
cx q[37], q[41];
t q[41];
cx q[37], q[55];
tdg q[55];
cx q[37], q[55];
t q[37];
t q[55];
h q[55];
t q[35];
cx q[87], q[99];
h q[81];
h q[51];
cx q[51], q[94];
tdg q[94];
cx q[54], q[94];
t q[94];
cx q[51], q[94];
tdg q[94];
cx q[54], q[94];
t q[94];
cx q[54], q[51];
tdg q[51];
cx q[54], q[51];
t q[54];
t q[51];
h q[51];
t q[99];
t q[48];
h q[19];
cx q[19], q[17];
tdg q[17];
cx q[58], q[17];
t q[17];
cx q[19], q[17];
tdg q[17];
cx q[58], q[17];
t q[17];
cx q[58], q[19];
tdg q[19];
cx q[58], q[19];
t q[58];
t q[19];
h q[19];
h q[2];
s q[37];
cx q[68], q[31];
t q[13];
h q[43];
h q[86];
cx q[86], q[49];
tdg q[49];
cx q[9], q[49];
t q[49];
cx q[86], q[49];
tdg q[49];
cx q[9], q[49];
t q[49];
cx q[9], q[86];
tdg q[86];
cx q[9], q[86];
t q[9];
t q[86];
h q[86];
cx q[89], q[50];
h q[24];
h q[56];
h q[62];
t q[93];
cx q[43], q[66];
s q[10];
h q[52];
s q[24];
h q[23];
cx q[23], q[5];
tdg q[5];
cx q[29], q[5];
t q[5];
cx q[23], q[5];
tdg q[5];
cx q[29], q[5];
t q[5];
cx q[29], q[23];
tdg q[23];
cx q[29], q[23];
t q[29];
t q[23];
h q[23];
h q[49];
cx q[49], q[77];
tdg q[77];
cx q[97], q[77];
t q[77];
cx q[49], q[77];
tdg q[77];
cx q[97], q[77];
t q[77];
cx q[97], q[49];
tdg q[49];
cx q[97], q[49];
t q[97];
t q[49];
h q[49];
h q[47];
cx q[47], q[84];
tdg q[84];
cx q[91], q[84];
t q[84];
cx q[47], q[84];
tdg q[84];
cx q[91], q[84];
t q[84];
cx q[91], q[47];
tdg q[47];
cx q[91], q[47];
t q[91];
t q[47];
h q[47];
cx q[99], q[82];
s q[64];
t q[70];
h q[1];
cx q[1], q[62];
tdg q[62];
cx q[36], q[62];
t q[62];
cx q[1], q[62];
tdg q[62];
cx q[36], q[62];
t q[62];
cx q[36], q[1];
tdg q[1];
cx q[36], q[1];
t q[36];
t q[1];
h q[1];
h q[22];
h q[9];
cx q[9], q[41];
tdg q[41];
cx q[90], q[41];
t q[41];
cx q[9], q[41];
tdg q[41];
cx q[90], q[41];
t q[41];
cx q[90], q[9];
tdg q[9];
cx q[90], q[9];
t q[90];
t q[9];
h q[9];
t q[28];
h q[24];
cx q[24], q[11];
tdg q[11];
cx q[55], q[11];
t q[11];
cx q[24], q[11];
tdg q[11];
cx q[55], q[11];
t q[11];
cx q[55], q[24];
tdg q[24];
cx q[55], q[24];
t q[55];
t q[24];
h q[24];
t q[30];
cx q[46], q[24];
s q[87];
h q[62];
cx q[62], q[95];
tdg q[95];
cx q[47], q[95];
t q[95];
cx q[62], q[95];
tdg q[95];
cx q[47], q[95];
t q[95];
cx q[47], q[62];
tdg q[62];
cx q[47], q[62];
t q[47];
t q[62];
h q[62];
s q[96];
cx q[92], q[36];
cx q[18], q[24];
s q[39];
s q[53];
h q[58];
cx q[58], q[23];
tdg q[23];
cx q[93], q[23];
t q[23];
cx q[58], q[23];
tdg q[23];
cx q[93], q[23];
t q[23];
cx q[93], q[58];
tdg q[58];
cx q[93], q[58];
t q[93];
t q[58];
h q[58];
cx q[80], q[26];
h q[3];
cx q[3], q[23];
tdg q[23];
cx q[19], q[23];
t q[23];
cx q[3], q[23];
tdg q[23];
cx q[19], q[23];
t q[23];
cx q[19], q[3];
tdg q[3];
cx q[19], q[3];
t q[19];
t q[3];
h q[3];
h q[43];
cx q[70], q[20];
h q[45];
cx q[0], q[60];
t q[12];
h q[22];
cx q[22], q[47];
tdg q[47];
cx q[95], q[47];
t q[47];
cx q[22], q[47];
tdg q[47];
cx q[95], q[47];
t q[47];
cx q[95], q[22];
tdg q[22];
cx q[95], q[22];
t q[95];
t q[22];
h q[22];
h q[30];
cx q[30], q[42];
tdg q[42];
cx q[39], q[42];
t q[42];
cx q[30], q[42];
tdg q[42];
cx q[39], q[42];
t q[42];
cx q[39], q[30];
tdg q[30];
cx q[39], q[30];
t q[39];
t q[30];
h q[30];
t q[29];
t q[68];
s q[10];
h q[75];
cx q[75], q[67];
tdg q[67];
cx q[43], q[67];
t q[67];
cx q[75], q[67];
tdg q[67];
cx q[43], q[67];
t q[67];
cx q[43], q[75];
tdg q[75];
cx q[43], q[75];
t q[43];
t q[75];
h q[75];
t q[64];
h q[42];
cx q[38], q[5];
h q[94];
s q[7];
t q[76];
cx q[79], q[83];
cx q[4], q[75];
t q[66];
h q[41];
h q[82];
cx q[82], q[5];
tdg q[5];
cx q[24], q[5];
t q[5];
cx q[82], q[5];
tdg q[5];
cx q[24], q[5];
t q[5];
cx q[24], q[82];
tdg q[82];
cx q[24], q[82];
t q[24];
t q[82];
h q[82];
h q[26];
h q[96];
h q[64];
cx q[64], q[93];
tdg q[93];
cx q[77], q[93];
t q[93];
cx q[64], q[93];
tdg q[93];
cx q[77], q[93];
t q[93];
cx q[77], q[64];
tdg q[64];
cx q[77], q[64];
t q[77];
t q[64];
h q[64];
s q[79];
cx q[42], q[85];
h q[94];
cx q[94], q[79];
tdg q[79];
cx q[41], q[79];
t q[79];
cx q[94], q[79];
tdg q[79];
cx q[41], q[79];
t q[79];
cx q[41], q[94];
tdg q[94];
cx q[41], q[94];
t q[41];
t q[94];
h q[94];
h q[18];
s q[8];
cx q[43], q[99];
s q[75];
cx q[18], q[10];
s q[91];
h q[71];
s q[35];
h q[49];
h q[21];
cx q[21], q[65];
tdg q[65];
cx q[60], q[65];
t q[65];
cx q[21], q[65];
tdg q[65];
cx q[60], q[65];
t q[65];
cx q[60], q[21];
tdg q[21];
cx q[60], q[21];
t q[60];
t q[21];
h q[21];
h q[67];
h q[20];
cx q[20], q[39];
tdg q[39];
cx q[43], q[39];
t q[39];
cx q[20], q[39];
tdg q[39];
cx q[43], q[39];
t q[39];
cx q[43], q[20];
tdg q[20];
cx q[43], q[20];
t q[43];
t q[20];
h q[20];
h q[28];
cx q[28], q[91];
tdg q[91];
cx q[44], q[91];
t q[91];
cx q[28], q[91];
tdg q[91];
cx q[44], q[91];
t q[91];
cx q[44], q[28];
tdg q[28];
cx q[44], q[28];
t q[44];
t q[28];
h q[28];
cx q[87], q[72];
t q[14];
cx q[30], q[72];
cx q[88], q[27];
h q[60];
t q[36];
h q[51];
cx q[51], q[4];
tdg q[4];
cx q[16], q[4];
t q[4];
cx q[51], q[4];
tdg q[4];
cx q[16], q[4];
t q[4];
cx q[16], q[51];
tdg q[51];
cx q[16], q[51];
t q[16];
t q[51];
h q[51];
h q[95];
cx q[95], q[14];
tdg q[14];
cx q[84], q[14];
t q[14];
cx q[95], q[14];
tdg q[14];
cx q[84], q[14];
t q[14];
cx q[84], q[95];
tdg q[95];
cx q[84], q[95];
t q[84];
t q[95];
h q[95];
cx q[21], q[17];
h q[99];
cx q[99], q[84];
tdg q[84];
cx q[51], q[84];
t q[84];
cx q[99], q[84];
tdg q[84];
cx q[51], q[84];
t q[84];
cx q[51], q[99];
tdg q[99];
cx q[51], q[99];
t q[51];
t q[99];
h q[99];
h q[63];
h q[52];
cx q[52], q[84];
tdg q[84];
cx q[92], q[84];
t q[84];
cx q[52], q[84];
tdg q[84];
cx q[92], q[84];
t q[84];
cx q[92], q[52];
tdg q[52];
cx q[92], q[52];
t q[92];
t q[52];
h q[52];
s q[79];
t q[86];
cx q[39], q[92];
cx q[15], q[24];
t q[66];
t q[37];
h q[85];
cx q[85], q[15];
tdg q[15];
cx q[11], q[15];
t q[15];
cx q[85], q[15];
tdg q[15];
cx q[11], q[15];
t q[15];
cx q[11], q[85];
tdg q[85];
cx q[11], q[85];
t q[11];
t q[85];
h q[85];
t q[18];
s q[80];
t q[34];
s q[59];
s q[64];
t q[91];
s q[25];
s q[85];
s q[94];
h q[13];
h q[4];
cx q[4], q[94];
tdg q[94];
cx q[29], q[94];
t q[94];
cx q[4], q[94];
tdg q[94];
cx q[29], q[94];
t q[94];
cx q[29], q[4];
tdg q[4];
cx q[29], q[4];
t q[29];
t q[4];
h q[4];
h q[93];
h q[32];
cx q[32], q[85];
tdg q[85];
cx q[56], q[85];
t q[85];
cx q[32], q[85];
tdg q[85];
cx q[56], q[85];
t q[85];
cx q[56], q[32];
tdg q[32];
cx q[56], q[32];
t q[56];
t q[32];
h q[32];
t q[72];
h q[28];
cx q[46], q[53];
s q[22];
cx q[28], q[79];
h q[0];
cx q[0], q[29];
tdg q[29];
cx q[31], q[29];
t q[29];
cx q[0], q[29];
tdg q[29];
cx q[31], q[29];
t q[29];
cx q[31], q[0];
tdg q[0];
cx q[31], q[0];
t q[31];
t q[0];
h q[0];
cx q[19], q[98];
h q[48];
cx q[48], q[91];
tdg q[91];
cx q[0], q[91];
t q[91];
cx q[48], q[91];
tdg q[91];
cx q[0], q[91];
t q[91];
cx q[0], q[48];
tdg q[48];
cx q[0], q[48];
t q[0];
t q[48];
h q[48];
cx q[65], q[97];
t q[26];
cx q[16], q[31];
cx q[50], q[46];
t q[68];
h q[18];
cx q[18], q[31];
tdg q[31];
cx q[42], q[31];
t q[31];
cx q[18], q[31];
tdg q[31];
cx q[42], q[31];
t q[31];
cx q[42], q[18];
tdg q[18];
cx q[42], q[18];
t q[42];
t q[18];
h q[18];
cx q[46], q[3];
h q[41];
cx q[41], q[64];
tdg q[64];
cx q[66], q[64];
t q[64];
cx q[41], q[64];
tdg q[64];
cx q[66], q[64];
t q[64];
cx q[66], q[41];
tdg q[41];
cx q[66], q[41];
t q[66];
t q[41];
h q[41];
h q[5];
cx q[5], q[31];
tdg q[31];
cx q[66], q[31];
t q[31];
cx q[5], q[31];
tdg q[31];
cx q[66], q[31];
t q[31];
cx q[66], q[5];
tdg q[5];
cx q[66], q[5];
t q[66];
t q[5];
h q[5];
h q[52];
cx q[52], q[54];
tdg q[54];
cx q[17], q[54];
t q[54];
cx q[52], q[54];
tdg q[54];
cx q[17], q[54];
t q[54];
cx q[17], q[52];
tdg q[52];
cx q[17], q[52];
t q[17];
t q[52];
h q[52];
h q[49];
cx q[49], q[28];
tdg q[28];
cx q[12], q[28];
t q[28];
cx q[49], q[28];
tdg q[28];
cx q[12], q[28];
t q[28];
cx q[12], q[49];
tdg q[49];
cx q[12], q[49];
t q[12];
t q[49];
h q[49];
s q[76];
s q[62];
h q[42];
cx q[42], q[58];
tdg q[58];
cx q[93], q[58];
t q[58];
cx q[42], q[58];
tdg q[58];
cx q[93], q[58];
t q[58];
cx q[93], q[42];
tdg q[42];
cx q[93], q[42];
t q[93];
t q[42];
h q[42];
h q[76];
h q[71];
cx q[71], q[35];
tdg q[35];
cx q[6], q[35];
t q[35];
cx q[71], q[35];
tdg q[35];
cx q[6], q[35];
t q[35];
cx q[6], q[71];
tdg q[71];
cx q[6], q[71];
t q[6];
t q[71];
h q[71];
cx q[19], q[24];
h q[73];
cx q[73], q[1];
tdg q[1];
cx q[12], q[1];
t q[1];
cx q[73], q[1];
tdg q[1];
cx q[12], q[1];
t q[1];
cx q[12], q[73];
tdg q[73];
cx q[12], q[73];
t q[12];
t q[73];
h q[73];
h q[22];
h q[28];
h q[6];
cx q[6], q[89];
tdg q[89];
cx q[72], q[89];
t q[89];
cx q[6], q[89];
tdg q[89];
cx q[72], q[89];
t q[89];
cx q[72], q[6];
tdg q[6];
cx q[72], q[6];
t q[72];
t q[6];
h q[6];
cx q[5], q[7];
s q[32];
t q[55];
h q[53];
s q[60];
cx q[65], q[87];
h q[16];
t q[36];
h q[57];
t q[27];
cx q[36], q[86];
cx q[32], q[2];
s q[9];
s q[71];
h q[72];
cx q[72], q[53];
tdg q[53];
cx q[15], q[53];
t q[53];
cx q[72], q[53];
tdg q[53];
cx q[15], q[53];
t q[53];
cx q[15], q[72];
tdg q[72];
cx q[15], q[72];
t q[15];
t q[72];
h q[72];
h q[24];
h q[91];
h q[90];
s q[69];
cx q[29], q[26];
h q[81];
h q[84];
cx q[21], q[44];
s q[71];
h q[5];
s q[42];
s q[50];
cx q[34], q[3];
t q[91];
t q[38];
h q[47];
h q[22];
h q[67];
cx q[67], q[11];
tdg q[11];
cx q[49], q[11];
t q[11];
cx q[67], q[11];
tdg q[11];
cx q[49], q[11];
t q[11];
cx q[49], q[67];
tdg q[67];
cx q[49], q[67];
t q[49];
t q[67];
h q[67];
cx q[44], q[93];
h q[95];
t q[65];
s q[61];
h q[9];
cx q[9], q[74];
tdg q[74];
cx q[51], q[74];
t q[74];
cx q[9], q[74];
tdg q[74];
cx q[51], q[74];
t q[74];
cx q[51], q[9];
tdg q[9];
cx q[51], q[9];
t q[51];
t q[9];
h q[9];
h q[86];
cx q[86], q[47];
tdg q[47];
cx q[19], q[47];
t q[47];
cx q[86], q[47];
tdg q[47];
cx q[19], q[47];
t q[47];
cx q[19], q[86];
tdg q[86];
cx q[19], q[86];
t q[19];
t q[86];
h q[86];
h q[22];
cx q[22], q[15];
tdg q[15];
cx q[4], q[15];
t q[15];
cx q[22], q[15];
tdg q[15];
cx q[4], q[15];
t q[15];
cx q[4], q[22];
tdg q[22];
cx q[4], q[22];
t q[4];
t q[22];
h q[22];
s q[4];
h q[95];
cx q[95], q[32];
tdg q[32];
cx q[83], q[32];
t q[32];
cx q[95], q[32];
tdg q[32];
cx q[83], q[32];
t q[32];
cx q[83], q[95];
tdg q[95];
cx q[83], q[95];
t q[83];
t q[95];
h q[95];
h q[35];
h q[46];
cx q[46], q[26];
tdg q[26];
cx q[57], q[26];
t q[26];
cx q[46], q[26];
tdg q[26];
cx q[57], q[26];
t q[26];
cx q[57], q[46];
tdg q[46];
cx q[57], q[46];
t q[57];
t q[46];
h q[46];
h q[3];
t q[82];
cx q[80], q[53];
h q[93];
h q[64];
cx q[64], q[72];
tdg q[72];
cx q[84], q[72];
t q[72];
cx q[64], q[72];
tdg q[72];
cx q[84], q[72];
t q[72];
cx q[84], q[64];
tdg q[64];
cx q[84], q[64];
t q[84];
t q[64];
h q[64];
t q[75];
t q[34];
h q[22];
cx q[22], q[55];
tdg q[55];
cx q[20], q[55];
t q[55];
cx q[22], q[55];
tdg q[55];
cx q[20], q[55];
t q[55];
cx q[20], q[22];
tdg q[22];
cx q[20], q[22];
t q[20];
t q[22];
h q[22];
h q[70];
cx q[70], q[89];
tdg q[89];
cx q[36], q[89];
t q[89];
cx q[70], q[89];
tdg q[89];
cx q[36], q[89];
t q[89];
cx q[36], q[70];
tdg q[70];
cx q[36], q[70];
t q[36];
t q[70];
h q[70];
h q[90];
cx q[90], q[64];
tdg q[64];
cx q[58], q[64];
t q[64];
cx q[90], q[64];
tdg q[64];
cx q[58], q[64];
t q[64];
cx q[58], q[90];
tdg q[90];
cx q[58], q[90];
t q[58];
t q[90];
h q[90];
cx q[6], q[79];
h q[55];
cx q[27], q[41];
s q[43];
h q[29];
h q[0];
cx q[0], q[58];
tdg q[58];
cx q[31], q[58];
t q[58];
cx q[0], q[58];
tdg q[58];
cx q[31], q[58];
t q[58];
cx q[31], q[0];
tdg q[0];
cx q[31], q[0];
t q[31];
t q[0];
h q[0];
t q[79];
h q[57];
h q[76];
cx q[76], q[79];
tdg q[79];
cx q[97], q[79];
t q[79];
cx q[76], q[79];
tdg q[79];
cx q[97], q[79];
t q[79];
cx q[97], q[76];
tdg q[76];
cx q[97], q[76];
t q[97];
t q[76];
h q[76];
h q[1];
cx q[44], q[32];
s q[80];
cx q[2], q[51];
h q[91];
cx q[91], q[51];
tdg q[51];
cx q[96], q[51];
t q[51];
cx q[91], q[51];
tdg q[51];
cx q[96], q[51];
t q[51];
cx q[96], q[91];
tdg q[91];
cx q[96], q[91];
t q[96];
t q[91];
h q[91];
h q[98];
cx q[98], q[23];
tdg q[23];
cx q[68], q[23];
t q[23];
cx q[98], q[23];
tdg q[23];
cx q[68], q[23];
t q[23];
cx q[68], q[98];
tdg q[98];
cx q[68], q[98];
t q[68];
t q[98];
h q[98];
h q[8];
cx q[8], q[16];
tdg q[16];
cx q[66], q[16];
t q[16];
cx q[8], q[16];
tdg q[16];
cx q[66], q[16];
t q[16];
cx q[66], q[8];
tdg q[8];
cx q[66], q[8];
t q[66];
t q[8];
h q[8];
t q[30];
s q[7];
s q[17];
h q[38];
cx q[38], q[44];
tdg q[44];
cx q[10], q[44];
t q[44];
cx q[38], q[44];
tdg q[44];
cx q[10], q[44];
t q[44];
cx q[10], q[38];
tdg q[38];
cx q[10], q[38];
t q[10];
t q[38];
h q[38];
h q[40];
cx q[40], q[47];
tdg q[47];
cx q[96], q[47];
t q[47];
cx q[40], q[47];
tdg q[47];
cx q[96], q[47];
t q[47];
cx q[96], q[40];
tdg q[40];
cx q[96], q[40];
t q[96];
t q[40];
h q[40];
cx q[40], q[92];
s q[57];
h q[24];
s q[46];
cx q[18], q[95];
s q[90];
h q[54];
cx q[54], q[79];
tdg q[79];
cx q[87], q[79];
t q[79];
cx q[54], q[79];
tdg q[79];
cx q[87], q[79];
t q[79];
cx q[87], q[54];
tdg q[54];
cx q[87], q[54];
t q[87];
t q[54];
h q[54];
cx q[96], q[26];
t q[19];
h q[36];
cx q[36], q[13];
tdg q[13];
cx q[39], q[13];
t q[13];
cx q[36], q[13];
tdg q[13];
cx q[39], q[13];
t q[13];
cx q[39], q[36];
tdg q[36];
cx q[39], q[36];
t q[39];
t q[36];
h q[36];
t q[76];
h q[32];
h q[42];
cx q[42], q[12];
tdg q[12];
cx q[97], q[12];
t q[12];
cx q[42], q[12];
tdg q[12];
cx q[97], q[12];
t q[12];
cx q[97], q[42];
tdg q[42];
cx q[97], q[42];
t q[97];
t q[42];
h q[42];
cx q[1], q[62];
cx q[81], q[21];
t q[41];
h q[81];
cx q[81], q[35];
tdg q[35];
cx q[43], q[35];
t q[35];
cx q[81], q[35];
tdg q[35];
cx q[43], q[35];
t q[35];
cx q[43], q[81];
tdg q[81];
cx q[43], q[81];
t q[43];
t q[81];
h q[81];
h q[26];
cx q[26], q[10];
tdg q[10];
cx q[12], q[10];
t q[10];
cx q[26], q[10];
tdg q[10];
cx q[12], q[10];
t q[10];
cx q[12], q[26];
tdg q[26];
cx q[12], q[26];
t q[12];
t q[26];
h q[26];
h q[83];
t q[57];
h q[95];
cx q[95], q[54];
tdg q[54];
cx q[83], q[54];
t q[54];
cx q[95], q[54];
tdg q[54];
cx q[83], q[54];
t q[54];
cx q[83], q[95];
tdg q[95];
cx q[83], q[95];
t q[83];
t q[95];
h q[95];
s q[71];
h q[36];
t q[49];
cx q[49], q[46];
cx q[16], q[19];
cx q[51], q[67];
h q[52];
t q[64];
h q[55];
cx q[94], q[87];
t q[47];
t q[80];
h q[89];
cx q[89], q[23];
tdg q[23];
cx q[43], q[23];
t q[23];
cx q[89], q[23];
tdg q[23];
cx q[43], q[23];
t q[23];
cx q[43], q[89];
tdg q[89];
cx q[43], q[89];
t q[43];
t q[89];
h q[89];
t q[87];
h q[92];
cx q[92], q[88];
tdg q[88];
cx q[32], q[88];
t q[88];
cx q[92], q[88];
tdg q[88];
cx q[32], q[88];
t q[88];
cx q[32], q[92];
tdg q[92];
cx q[32], q[92];
t q[32];
t q[92];
h q[92];
t q[55];
t q[23];
s q[58];
t q[50];
cx q[55], q[98];
h q[7];
h q[36];
h q[26];
cx q[26], q[81];
tdg q[81];
cx q[76], q[81];
t q[81];
cx q[26], q[81];
tdg q[81];
cx q[76], q[81];
t q[81];
cx q[76], q[26];
tdg q[26];
cx q[76], q[26];
t q[76];
t q[26];
h q[26];
cx q[62], q[26];
h q[30];
cx q[74], q[19];
h q[66];
h q[68];
cx q[68], q[67];
tdg q[67];
cx q[74], q[67];
t q[67];
cx q[68], q[67];
tdg q[67];
cx q[74], q[67];
t q[67];
cx q[74], q[68];
tdg q[68];
cx q[74], q[68];
t q[74];
t q[68];
h q[68];
cx q[81], q[89];
t q[49];
h q[0];
cx q[0], q[25];
tdg q[25];
cx q[4], q[25];
t q[25];
cx q[0], q[25];
tdg q[25];
cx q[4], q[25];
t q[25];
cx q[4], q[0];
tdg q[0];
cx q[4], q[0];
t q[4];
t q[0];
h q[0];
s q[41];
cx q[10], q[9];
t q[19];
h q[50];
h q[58];
h q[66];
cx q[66], q[73];
tdg q[73];
cx q[90], q[73];
t q[73];
cx q[66], q[73];
tdg q[73];
cx q[90], q[73];
t q[73];
cx q[90], q[66];
tdg q[66];
cx q[90], q[66];
t q[90];
t q[66];
h q[66];
cx q[14], q[22];
s q[83];
h q[53];
cx q[53], q[0];
tdg q[0];
cx q[16], q[0];
t q[0];
cx q[53], q[0];
tdg q[0];
cx q[16], q[0];
t q[0];
cx q[16], q[53];
tdg q[53];
cx q[16], q[53];
t q[16];
t q[53];
h q[53];
s q[27];
cx q[2], q[3];
s q[3];
cx q[2], q[3];
cx q[5], q[4];
s q[4];
cx q[5], q[4];
s q[4];
t q[6];
t q[7];
s q[6];
cx q[6], q[7];
cx q[9], q[8];
y q[10];
z q[13];
x q[15];
s q[19];
t q[19];
cx q[18], q[19];
t q[21];
s q[20];
cx q[21], q[20];
cx q[23], q[22];
t q[22];
s q[22];
cx q[23], q[22];
t q[24];
s q[25];
cx q[24], q[25];
t q[26];
t q[27];
cx q[27], q[26];
cx q[29], q[28];
s q[28];
cx q[29], q[28];
cx q[30], q[31];
z q[32];
x q[34];
cx q[36], q[35];
s q[35];
cx q[36], q[35];
t q[35];
cx q[42], q[43];
z q[44];
cx q[46], q[45];
s q[48];
cx q[47], q[48];
s q[48];
cx q[47], q[48];
t q[93];
t q[57];
t q[91];
s q[54];
t q[50];
h q[75];
h q[72];
h q[51];
cx q[51], q[81];
tdg q[81];
cx q[89], q[81];
t q[81];
cx q[51], q[81];
tdg q[81];
cx q[89], q[81];
t q[81];
cx q[89], q[51];
tdg q[51];
cx q[89], q[51];
t q[89];
t q[51];
h q[51];
cx q[78], q[51];
h q[82];
h q[71];
cx q[71], q[91];
tdg q[91];
cx q[89], q[91];
t q[91];
cx q[71], q[91];
tdg q[91];
cx q[89], q[91];
t q[91];
cx q[89], q[71];
tdg q[71];
cx q[89], q[71];
t q[89];
t q[71];
h q[71];
cx q[67], q[91];
s q[78];
cx q[54], q[53];
h q[96];
h q[75];
cx q[75], q[55];
tdg q[55];
cx q[78], q[55];
t q[55];
cx q[75], q[55];
tdg q[55];
cx q[78], q[55];
t q[55];
cx q[78], q[75];
tdg q[75];
cx q[78], q[75];
t q[78];
t q[75];
h q[75];
h q[91];
cx q[91], q[78];
tdg q[78];
cx q[61], q[78];
t q[78];
cx q[91], q[78];
tdg q[78];
cx q[61], q[78];
t q[78];
cx q[61], q[91];
tdg q[91];
cx q[61], q[91];
t q[61];
t q[91];
h q[91];
h q[52];
cx q[52], q[72];
tdg q[72];
cx q[97], q[72];
t q[72];
cx q[52], q[72];
tdg q[72];
cx q[97], q[72];
t q[72];
cx q[97], q[52];
tdg q[52];
cx q[97], q[52];
t q[97];
t q[52];
h q[52];
h q[51];
cx q[51], q[52];
tdg q[52];
cx q[94], q[52];
t q[52];
cx q[51], q[52];
tdg q[52];
cx q[94], q[52];
t q[52];
cx q[94], q[51];
tdg q[51];
cx q[94], q[51];
t q[94];
t q[51];
h q[51];
h q[51];
h q[70];
h q[98];
h q[73];
cx q[73], q[79];
tdg q[79];
cx q[63], q[79];
t q[79];
cx q[73], q[79];
tdg q[79];
cx q[63], q[79];
t q[79];
cx q[63], q[73];
tdg q[73];
cx q[63], q[73];
t q[63];
t q[73];
h q[73];
h q[53];
cx q[53], q[69];
tdg q[69];
cx q[70], q[69];
t q[69];
cx q[53], q[69];
tdg q[69];
cx q[70], q[69];
t q[69];
cx q[70], q[53];
tdg q[53];
cx q[70], q[53];
t q[70];
t q[53];
h q[53];
cx q[64], q[82];
s q[57];
cx q[66], q[81];
s q[83];
h q[93];
s q[85];
s q[78];
h q[91];
s q[67];
h q[74];
h q[66];
h q[82];
cx q[82], q[80];
tdg q[80];
cx q[92], q[80];
t q[80];
cx q[82], q[80];
tdg q[80];
cx q[92], q[80];
t q[80];
cx q[92], q[82];
tdg q[82];
cx q[92], q[82];
t q[92];
t q[82];
h q[82];
t q[70];
h q[77];
cx q[77], q[98];
tdg q[98];
cx q[74], q[98];
t q[98];
cx q[77], q[98];
tdg q[98];
cx q[74], q[98];
t q[98];
cx q[74], q[77];
tdg q[77];
cx q[74], q[77];
t q[74];
t q[77];
h q[77];
s q[80];
h q[50];
cx q[50], q[95];
tdg q[95];
cx q[80], q[95];
t q[95];
cx q[50], q[95];
tdg q[95];
cx q[80], q[95];
t q[95];
cx q[80], q[50];
tdg q[50];
cx q[80], q[50];
t q[80];
t q[50];
h q[50];
h q[56];
cx q[69], q[78];
cx q[62], q[63];
h q[84];
cx q[84], q[60];
tdg q[60];
cx q[66], q[60];
t q[60];
cx q[84], q[60];
tdg q[60];
cx q[66], q[60];
t q[60];
cx q[66], q[84];
tdg q[84];
cx q[66], q[84];
t q[66];
t q[84];
h q[84];
s q[53];
h q[58];
cx q[59], q[77];
h q[88];
t q[82];
h q[66];
cx q[66], q[55];
tdg q[55];
cx q[69], q[55];
t q[55];
cx q[66], q[55];
tdg q[55];
cx q[69], q[55];
t q[55];
cx q[69], q[66];
tdg q[66];
cx q[69], q[66];
t q[69];
t q[66];
h q[66];
