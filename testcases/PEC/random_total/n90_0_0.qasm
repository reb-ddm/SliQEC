OPENQASM 2.0;
include "qelib1.inc";
qreg q[90];
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
t q[52];
ccx q[48], q[69], q[87];
h q[33];
ccx q[65], q[41], q[18];
cx q[84], q[3];
s q[48];
cx q[45], q[83];
h q[52];
ccx q[21], q[10], q[75];
cx q[11], q[86];
ccx q[65], q[44], q[76];
ccx q[4], q[84], q[22];
h q[24];
cx q[63], q[73];
h q[44];
ccx q[58], q[65], q[56];
cx q[84], q[85];
h q[24];
ccx q[25], q[61], q[86];
t q[61];
t q[85];
ccx q[5], q[61], q[33];
cx q[11], q[79];
t q[8];
s q[28];
cx q[88], q[52];
ccx q[72], q[48], q[13];
cx q[58], q[32];
s q[12];
h q[80];
s q[66];
h q[3];
ccx q[48], q[56], q[5];
h q[24];
s q[35];
cx q[19], q[23];
h q[75];
ccx q[73], q[61], q[62];
ccx q[13], q[2], q[24];
s q[87];
h q[48];
ccx q[53], q[29], q[16];
s q[32];
t q[24];
s q[46];
ccx q[75], q[39], q[36];
s q[75];
t q[53];
cx q[83], q[26];
h q[89];
cx q[28], q[54];
t q[51];
s q[15];
cx q[15], q[43];
s q[83];
s q[79];
h q[66];
cx q[45], q[18];
ccx q[35], q[53], q[66];
s q[15];
ccx q[0], q[81], q[71];
t q[46];
t q[52];
ccx q[25], q[84], q[13];
h q[38];
ccx q[33], q[87], q[63];
h q[17];
ccx q[72], q[86], q[6];
cx q[56], q[26];
s q[54];
s q[13];
cx q[53], q[54];
h q[14];
t q[73];
h q[83];
s q[0];
s q[19];
cx q[88], q[78];
t q[12];
cx q[34], q[3];
cx q[19], q[52];
t q[12];
h q[33];
ccx q[63], q[43], q[45];
cx q[55], q[64];
cx q[86], q[4];
t q[78];
cx q[72], q[40];
cx q[39], q[82];
t q[27];
h q[1];
ccx q[4], q[67], q[84];
t q[33];
h q[15];
ccx q[75], q[5], q[28];
cx q[49], q[22];
cx q[1], q[81];
t q[11];
t q[85];
ccx q[50], q[73], q[34];
ccx q[35], q[65], q[60];
ccx q[54], q[43], q[24];
s q[85];
t q[83];
cx q[15], q[29];
cx q[66], q[18];
ccx q[17], q[85], q[27];
ccx q[21], q[69], q[63];
h q[41];
h q[83];
cx q[59], q[35];
t q[5];
cx q[78], q[12];
s q[45];
h q[43];
ccx q[21], q[14], q[70];
h q[13];
ccx q[4], q[87], q[42];
s q[34];
cx q[85], q[12];
s q[33];
s q[6];
t q[9];
h q[28];
s q[80];
h q[35];
cx q[84], q[52];
ccx q[41], q[72], q[53];
t q[9];
s q[38];
s q[5];
cx q[13], q[6];
h q[19];
s q[36];
s q[28];
h q[61];
s q[84];
s q[63];
t q[34];
cx q[74], q[54];
s q[70];
ccx q[9], q[73], q[4];
h q[89];
ccx q[76], q[6], q[21];
t q[23];
t q[36];
cx q[50], q[35];
t q[13];
ccx q[61], q[79], q[85];
ccx q[55], q[24], q[15];
h q[30];
ccx q[41], q[4], q[89];
s q[22];
cx q[82], q[77];
h q[73];
t q[79];
cx q[25], q[37];
ccx q[52], q[72], q[62];
s q[41];
cx q[53], q[44];
h q[76];
cx q[14], q[84];
t q[46];
h q[5];
t q[32];
h q[37];
h q[59];
s q[64];
cx q[84], q[83];
t q[32];
cx q[30], q[35];
h q[49];
h q[45];
h q[48];
t q[74];
cx q[42], q[80];
cx q[19], q[63];
s q[1];
h q[51];
cx q[23], q[30];
cx q[16], q[83];
h q[15];
h q[42];
ccx q[41], q[56], q[63];
t q[16];
cx q[22], q[8];
cx q[0], q[24];
h q[18];
s q[29];
h q[47];
s q[28];
t q[36];
ccx q[16], q[67], q[44];
cx q[29], q[43];
h q[89];
s q[43];
ccx q[38], q[82], q[41];
s q[17];
ccx q[27], q[7], q[4];
s q[85];
cx q[83], q[29];
s q[36];
s q[7];
ccx q[60], q[82], q[31];
s q[87];
ccx q[12], q[70], q[75];
s q[41];
ccx q[54], q[55], q[63];
t q[10];
ccx q[19], q[15], q[11];
s q[78];
s q[61];
cx q[39], q[7];
ccx q[9], q[73], q[87];
h q[51];
s q[36];
h q[86];
s q[44];
h q[57];
t q[46];
t q[66];
s q[41];
s q[40];
h q[48];
s q[22];
ccx q[34], q[62], q[29];
h q[45];
h q[86];
h q[10];
h q[31];
t q[42];
t q[49];
t q[69];
ccx q[36], q[19], q[35];
h q[8];
cx q[13], q[74];
s q[60];
s q[68];
t q[60];
h q[55];
h q[67];
t q[21];
t q[21];
s q[31];
ccx q[83], q[46], q[19];
t q[33];
cx q[74], q[35];
ccx q[50], q[63], q[30];
h q[51];
ccx q[67], q[26], q[10];
s q[16];
ccx q[31], q[4], q[86];
t q[48];
cx q[78], q[74];
h q[28];
t q[69];
t q[26];
s q[43];
ccx q[18], q[43], q[41];
ccx q[56], q[17], q[72];
t q[55];
cx q[50], q[34];
t q[75];
h q[22];
cx q[27], q[57];
t q[83];
ccx q[41], q[56], q[47];
ccx q[35], q[54], q[30];
h q[81];
h q[47];
