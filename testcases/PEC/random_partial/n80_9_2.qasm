OPENQASM 2.0;
include "qelib1.inc";
qreg q[92];
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
t q[35];
s q[53];
t q[49];
t q[41];
s q[7];
cx q[41], q[55];
cx q[36], q[14];
s q[65];
t q[35];
cx q[1], q[18];
cx q[63], q[24];
cx q[12], q[48];
ccx q[48], q[21], q[29];
s q[16];
t q[22];
cx q[32], q[66];
s q[47];
cx q[51], q[11];
h q[26];
s q[69];
h q[78];
cx q[2], q[15];
s q[61];
ccx q[69], q[12], q[50];
s q[73];
t q[4];
s q[53];
s q[46];
t q[4];
s q[47];
s q[45];
t q[1];
s q[43];
s q[64];
s q[76];
t q[42];
t q[21];
cx q[73], q[52];
s q[63];
s q[7];
s q[55];
cx q[24], q[65];
t q[68];
h q[60];
t q[47];
h q[31];
t q[78];
ccx q[79], q[30], q[32];
ccx q[13], q[26], q[8];
h q[71];
t q[21];
s q[6];
t q[40];
h q[57];
h q[74];
ccx q[30], q[31], q[15];
t q[7];
h q[53];
s q[56];
ccx q[60], q[36], q[56];
cx q[5], q[2];
s q[18];
s q[33];
ccx q[11], q[70], q[64];
t q[68];
ccx q[34], q[20], q[15];
t q[35];
ccx q[26], q[31], q[63];
cx q[42], q[29];
ccx q[59], q[67], q[35];
s q[30];
ccx q[68], q[21], q[77];
h q[14];
ccx q[71], q[58], q[67];
s q[13];
ccx q[68], q[72], q[13];
s q[38];
cx q[22], q[42];
s q[2];
cx q[58], q[28];
cx q[22], q[50];
h q[37];
s q[32];
s q[59];
t q[22];
ccx q[56], q[47], q[18];
s q[61];
ccx q[15], q[47], q[10];
cx q[59], q[73];
h q[29];
ccx q[74], q[0], q[22];
h q[54];
ccx q[70], q[25], q[40];
ccx q[53], q[22], q[50];
s q[67];
cx q[73], q[39];
cx q[49], q[12];
cx q[34], q[67];
cx q[54], q[5];
ccx q[77], q[23], q[71];
t q[59];
t q[48];
s q[43];
t q[2];
ccx q[26], q[10], q[23];
t q[58];
t q[52];
s q[73];
cx q[52], q[41];
ccx q[41], q[31], q[8];
ccx q[53], q[74], q[19];
h q[6];
ccx q[67], q[37], q[4];
h q[16];
t q[31];
h q[75];
h q[63];
t q[36];
cx q[21], q[25];
s q[3];
h q[62];
s q[19];
ccx q[46], q[58], q[31];
ccx q[3], q[65], q[23];
cx q[63], q[72];
ccx q[52], q[42], q[5];
ccx q[32], q[15], q[60];
ccx q[46], q[68], q[50];
ccx q[33], q[19], q[4];
t q[7];
h q[31];
s q[56];
ccx q[53], q[52], q[58];
s q[43];
h q[41];
h q[5];
h q[68];
ccx q[77], q[49], q[36];
s q[23];
cx q[75], q[76];
s q[13];
t q[73];
s q[10];
ccx q[53], q[15], q[61];
h q[12];
h q[49];
ccx q[62], q[8], q[27];
h q[21];
h q[18];
ccx q[39], q[10], q[25];
s q[51];
h q[20];
t q[36];
h q[45];
t q[52];
h q[62];
t q[31];
cx q[55], q[19];
ccx q[25], q[0], q[27];
s q[26];
s q[24];
s q[21];
cx q[44], q[8];
cx q[69], q[29];
t q[24];
ccx q[54], q[72], q[23];
ccx q[1], q[78], q[22];
t q[61];
t q[37];
ccx q[7], q[59], q[29];
s q[23];
cx q[33], q[4];
ccx q[24], q[44], q[55];
t q[44];
s q[41];
s q[52];
cx q[70], q[64];
h q[57];
ccx q[67], q[71], q[44];
s q[7];
h q[44];
cx q[10], q[8];
ccx q[33], q[17], q[26];
ccx q[51], q[37], q[55];
h q[45];
cx q[13], q[56];
h q[67];
cx q[63], q[45];
ccx q[0], q[45], q[37];
ccx q[3], q[68], q[75];
t q[55];
t q[13];
s q[47];
s q[29];
cx q[71], q[9];
cx q[27], q[78];
s q[14];
ccx q[46], q[0], q[39];
s q[78];
ccx q[45], q[13], q[24];
ccx q[31], q[71], q[52];
s q[39];
cx q[7], q[58];
t q[42];
t q[75];
cx q[44], q[34];
ccx q[53], q[75], q[4];
ccx q[4], q[45], q[57];
s q[20];
ccx q[77], q[78], q[63];
h q[63];
h q[28];
h q[6];
h q[27];
h q[6];
t q[39];
t q[34];
ccx q[36], q[63], q[35];
cx q[0], q[14];
t q[6];
h q[9];
t q[3];
h q[71];
ccx q[8], q[37], q[77];
s q[53];
s q[7];
cx q[59], q[70];
cx q[3], q[6];
ccx q[50], q[44], q[13];
t q[13];
cx q[57], q[39];
s q[15];
s q[63];
h q[74];
cx q[60], q[64];
h q[20];
t q[79];
s q[70];
cx q[43], q[23];
h q[60];
cx q[0], q[1];
sdg q[1];
cx q[0], q[1];
tdg q[0];
x q[2];
sdg q[3];
cx q[4], q[3];
sdg q[4];
tdg q[3];
z q[5];
y q[6];
sdg q[7];
z q[9];
sdg q[11];
cx q[11], q[10];
sdg q[10];
tdg q[12];
sdg q[14];
cx q[14], q[15];
tdg q[15];
sdg q[14];
cx q[14], q[15];
z q[16];
sdg q[17];
cx q[17], q[18];
x q[19];
cx q[20], q[21];
sdg q[20];
tdg q[23];
tdg q[22];
cx q[23], q[22];
tdg q[22];
y q[24];
sdg q[26];
cx q[27], q[26];
sdg q[27];
sdg q[26];
cx q[29], q[28];
x q[28];
cx q[30], q[31];
sdg q[32];
sdg q[34];
y q[36];
cx q[37], q[38];
tdg q[37];
tdg q[38];
cx q[37], q[38];
tdg q[37];
x q[39];
s q[66];
ccx q[51], q[62], q[42];
ccx q[74], q[79], q[66];
s q[74];
cx q[63], q[73];
ccx q[47], q[66], q[76];
t q[57];
cx q[42], q[79];
h q[75];
ccx q[45], q[56], q[44];
h q[78];
s q[43];
s q[73];
t q[76];
cx q[69], q[45];
ccx q[77], q[59], q[51];
s q[50];
cx q[43], q[69];
h q[46];
t q[62];
ccx q[71], q[78], q[67];
ccx q[48], q[61], q[52];
ccx q[51], q[55], q[50];
h q[42];
s q[40];
s q[74];
h q[75];
t q[65];
ccx q[51], q[55], q[49];
t q[53];
t q[73];
h q[68];
t q[68];
t q[43];
h q[77];
h q[77];
s q[42];
h q[70];
cx q[72], q[71];
cx q[74], q[67];
cx q[80], q[9];
cx q[81], q[27];
cx q[82], q[51];
cx q[83], q[0];
cx q[84], q[29];
cx q[85], q[58];
cx q[86], q[21];
cx q[87], q[34];
cx q[88], q[44];
cx q[89], q[46];
cx q[90], q[8];
cx q[91], q[10];
