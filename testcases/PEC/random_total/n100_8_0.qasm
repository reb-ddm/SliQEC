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
cx q[86], q[79];
t q[25];
t q[19];
cx q[12], q[98];
s q[28];
ccx q[0], q[19], q[44];
ccx q[18], q[44], q[52];
h q[14];
h q[93];
cx q[58], q[43];
cx q[25], q[60];
ccx q[69], q[30], q[58];
h q[92];
h q[81];
s q[3];
cx q[34], q[22];
ccx q[92], q[4], q[0];
ccx q[65], q[42], q[9];
cx q[97], q[46];
ccx q[77], q[49], q[62];
t q[33];
t q[88];
h q[22];
s q[1];
s q[16];
h q[74];
h q[68];
cx q[56], q[11];
s q[89];
cx q[20], q[86];
ccx q[44], q[28], q[1];
h q[79];
cx q[87], q[77];
t q[0];
t q[77];
h q[27];
h q[30];
cx q[53], q[32];
cx q[34], q[50];
h q[77];
t q[78];
t q[79];
t q[54];
cx q[27], q[90];
ccx q[55], q[59], q[16];
s q[27];
h q[40];
t q[98];
s q[69];
t q[4];
cx q[23], q[3];
ccx q[6], q[56], q[36];
h q[39];
ccx q[95], q[32], q[28];
s q[14];
h q[85];
s q[86];
cx q[8], q[35];
ccx q[16], q[45], q[51];
s q[11];
s q[36];
s q[94];
s q[83];
h q[42];
h q[86];
cx q[5], q[28];
cx q[12], q[78];
t q[46];
cx q[58], q[95];
cx q[9], q[74];
t q[10];
s q[60];
ccx q[82], q[50], q[34];
t q[14];
s q[32];
ccx q[88], q[90], q[40];
h q[86];
cx q[14], q[89];
t q[36];
ccx q[53], q[56], q[3];
ccx q[20], q[71], q[19];
ccx q[59], q[37], q[43];
cx q[11], q[85];
t q[22];
s q[32];
ccx q[97], q[74], q[30];
t q[11];
h q[66];
s q[75];
t q[78];
h q[39];
ccx q[23], q[32], q[46];
h q[0];
ccx q[15], q[72], q[70];
cx q[47], q[54];
t q[69];
cx q[21], q[2];
ccx q[32], q[90], q[51];
h q[32];
t q[26];
s q[88];
s q[21];
cx q[38], q[36];
t q[14];
h q[47];
t q[8];
cx q[0], q[76];
t q[4];
h q[41];
ccx q[96], q[6], q[83];
ccx q[14], q[2], q[41];
cx q[70], q[18];
cx q[92], q[7];
h q[91];
h q[56];
t q[91];
h q[58];
h q[11];
cx q[67], q[5];
cx q[85], q[94];
t q[3];
t q[44];
t q[34];
ccx q[33], q[50], q[25];
cx q[84], q[89];
cx q[49], q[15];
ccx q[17], q[66], q[74];
s q[44];
ccx q[67], q[5], q[2];
s q[81];
h q[75];
s q[51];
s q[42];
ccx q[27], q[6], q[88];
ccx q[92], q[11], q[34];
s q[49];
cx q[64], q[26];
cx q[64], q[88];
s q[64];
h q[41];
ccx q[51], q[40], q[52];
t q[93];
h q[58];
ccx q[51], q[62], q[6];
h q[94];
cx q[14], q[36];
s q[97];
cx q[58], q[59];
t q[54];
t q[47];
h q[94];
cx q[27], q[29];
cx q[89], q[83];
t q[51];
h q[11];
h q[19];
ccx q[7], q[14], q[37];
s q[79];
cx q[52], q[68];
h q[95];
ccx q[49], q[26], q[59];
h q[81];
h q[21];
t q[32];
ccx q[87], q[33], q[93];
ccx q[37], q[72], q[73];
cx q[5], q[14];
cx q[95], q[27];
ccx q[82], q[48], q[60];
ccx q[80], q[36], q[45];
t q[15];
s q[11];
s q[55];
t q[77];
h q[22];
s q[65];
ccx q[13], q[59], q[57];
cx q[0], q[4];
t q[36];
cx q[44], q[84];
h q[29];
s q[24];
ccx q[70], q[26], q[31];
s q[54];
cx q[33], q[41];
s q[15];
t q[10];
t q[88];
ccx q[74], q[70], q[95];
ccx q[43], q[0], q[48];
h q[22];
ccx q[26], q[25], q[47];
ccx q[58], q[96], q[78];
s q[2];
s q[81];
h q[20];
h q[64];
ccx q[57], q[39], q[88];
cx q[86], q[40];
cx q[12], q[18];
h q[54];
s q[35];
s q[3];
t q[70];
ccx q[91], q[81], q[93];
t q[78];
h q[9];
ccx q[13], q[74], q[45];
h q[41];
cx q[66], q[91];
ccx q[82], q[19], q[87];
ccx q[43], q[75], q[47];
t q[73];
t q[49];
s q[76];
h q[53];
t q[97];
t q[42];
ccx q[77], q[8], q[17];
cx q[89], q[7];
t q[78];
h q[44];
cx q[25], q[8];
cx q[6], q[94];
cx q[97], q[87];
t q[17];
cx q[22], q[99];
t q[48];
cx q[69], q[46];
cx q[10], q[45];
s q[94];
t q[85];
h q[32];
t q[21];
s q[69];
h q[60];
cx q[30], q[97];
cx q[38], q[55];
t q[74];
ccx q[12], q[70], q[31];
cx q[95], q[60];
ccx q[15], q[80], q[83];
s q[91];
s q[27];
cx q[79], q[46];
s q[96];
t q[62];
t q[40];
ccx q[33], q[23], q[47];
h q[59];
s q[34];
s q[33];
ccx q[37], q[11], q[79];
ccx q[53], q[1], q[29];
t q[97];
h q[40];
h q[33];
cx q[28], q[39];
h q[10];
t q[46];
h q[44];
h q[34];
h q[22];
s q[8];
t q[67];
t q[97];
t q[5];
t q[76];
t q[91];
cx q[85], q[9];
t q[26];
cx q[29], q[7];
s q[17];
ccx q[82], q[6], q[53];
ccx q[51], q[0], q[85];
ccx q[89], q[17], q[25];
h q[95];
ccx q[48], q[67], q[16];
h q[52];
cx q[5], q[70];
h q[70];
ccx q[59], q[22], q[15];
h q[50];
s q[16];
h q[38];
t q[81];
ccx q[37], q[18], q[85];
cx q[13], q[51];
t q[26];
cx q[56], q[81];
h q[84];
ccx q[10], q[86], q[70];
cx q[91], q[23];
h q[23];
s q[6];
t q[55];
t q[78];
ccx q[28], q[61], q[80];
cx q[50], q[21];
ccx q[50], q[81], q[26];
