OPENQASM 2.0;
include "qelib1.inc";
qreg q[106];
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
s q[55];
cx q[48], q[16];
cx q[20], q[81];
h q[81];
cx q[49], q[61];
t q[70];
h q[48];
cx q[82], q[55];
t q[16];
h q[52];
ccx q[83], q[4], q[63];
t q[3];
ccx q[23], q[51], q[49];
cx q[65], q[27];
t q[46];
t q[9];
h q[40];
s q[87];
s q[57];
cx q[46], q[13];
h q[66];
cx q[8], q[10];
s q[2];
s q[34];
t q[50];
cx q[53], q[49];
ccx q[54], q[78], q[60];
s q[29];
s q[10];
h q[58];
t q[22];
cx q[42], q[4];
t q[0];
ccx q[11], q[22], q[39];
h q[54];
t q[34];
ccx q[26], q[62], q[56];
cx q[24], q[6];
h q[89];
h q[60];
cx q[16], q[10];
cx q[2], q[5];
s q[23];
cx q[58], q[20];
t q[14];
s q[78];
h q[56];
cx q[78], q[77];
s q[80];
h q[36];
t q[51];
h q[36];
s q[59];
ccx q[35], q[79], q[18];
cx q[35], q[49];
cx q[84], q[50];
h q[24];
t q[26];
h q[61];
cx q[55], q[83];
t q[54];
cx q[65], q[69];
cx q[62], q[88];
ccx q[12], q[84], q[24];
cx q[45], q[5];
h q[0];
cx q[8], q[25];
t q[42];
h q[85];
h q[49];
s q[36];
t q[22];
ccx q[15], q[5], q[72];
ccx q[48], q[27], q[23];
cx q[64], q[26];
cx q[40], q[37];
s q[88];
h q[33];
t q[27];
cx q[27], q[5];
t q[62];
ccx q[60], q[83], q[18];
cx q[5], q[14];
t q[42];
cx q[42], q[50];
s q[8];
s q[67];
t q[12];
t q[29];
t q[44];
t q[32];
cx q[20], q[65];
cx q[55], q[79];
t q[27];
t q[11];
h q[84];
s q[31];
ccx q[75], q[89], q[32];
ccx q[44], q[27], q[57];
t q[4];
s q[50];
s q[71];
t q[2];
t q[58];
s q[54];
h q[10];
t q[47];
cx q[52], q[66];
ccx q[67], q[35], q[7];
t q[11];
ccx q[85], q[20], q[44];
t q[8];
cx q[45], q[12];
s q[19];
cx q[62], q[35];
t q[72];
h q[29];
ccx q[78], q[47], q[17];
h q[87];
t q[60];
t q[70];
t q[72];
h q[51];
s q[52];
h q[17];
ccx q[27], q[50], q[34];
h q[5];
t q[62];
t q[73];
s q[1];
s q[64];
ccx q[40], q[46], q[0];
s q[19];
h q[74];
cx q[14], q[17];
h q[14];
t q[39];
ccx q[20], q[52], q[61];
s q[11];
h q[76];
h q[11];
s q[16];
t q[31];
h q[73];
t q[83];
ccx q[55], q[23], q[11];
s q[56];
cx q[14], q[31];
h q[37];
cx q[30], q[77];
t q[0];
s q[35];
ccx q[24], q[72], q[69];
ccx q[80], q[12], q[46];
s q[60];
t q[34];
cx q[24], q[46];
s q[25];
cx q[14], q[81];
ccx q[77], q[82], q[84];
t q[38];
t q[28];
t q[50];
t q[73];
ccx q[32], q[60], q[66];
t q[46];
h q[40];
h q[59];
s q[5];
t q[29];
ccx q[6], q[78], q[0];
s q[81];
s q[62];
cx q[46], q[83];
t q[16];
s q[52];
t q[33];
cx q[24], q[61];
h q[81];
s q[10];
h q[14];
ccx q[56], q[40], q[78];
h q[59];
cx q[59], q[46];
ccx q[9], q[49], q[0];
cx q[71], q[40];
ccx q[81], q[32], q[58];
ccx q[62], q[0], q[25];
ccx q[20], q[32], q[23];
s q[35];
s q[10];
h q[43];
h q[59];
h q[21];
ccx q[56], q[40], q[36];
t q[9];
cx q[41], q[33];
ccx q[36], q[79], q[14];
h q[54];
t q[4];
s q[47];
ccx q[23], q[43], q[6];
h q[67];
s q[21];
s q[28];
h q[79];
t q[38];
cx q[1], q[19];
ccx q[75], q[74], q[62];
cx q[43], q[17];
cx q[51], q[55];
s q[0];
s q[16];
t q[84];
cx q[29], q[69];
ccx q[30], q[39], q[71];
s q[20];
cx q[71], q[73];
cx q[77], q[64];
t q[40];
t q[75];
t q[63];
ccx q[3], q[53], q[4];
cx q[73], q[43];
ccx q[37], q[4], q[17];
t q[76];
cx q[51], q[58];
h q[72];
s q[3];
h q[16];
cx q[70], q[1];
s q[34];
cx q[60], q[45];
ccx q[36], q[68], q[11];
t q[58];
h q[33];
t q[75];
t q[37];
cx q[86], q[46];
ccx q[27], q[67], q[7];
cx q[5], q[36];
s q[0];
t q[65];
h q[37];
ccx q[74], q[18], q[42];
cx q[79], q[84];
h q[38];
cx q[17], q[41];
h q[20];
h q[47];
h q[40];
t q[20];
h q[23];
s q[40];
cx q[40], q[68];
h q[57];
cx q[12], q[85];
cx q[44], q[46];
cx q[48], q[6];
cx q[12], q[83];
h q[47];
s q[37];
cx q[82], q[40];
s q[47];
ccx q[42], q[24], q[76];
cx q[72], q[86];
ccx q[21], q[87], q[35];
t q[56];
h q[67];
t q[35];
cx q[0], q[1];
tdg q[1];
sdg q[1];
cx q[0], q[1];
sdg q[0];
tdg q[3];
tdg q[2];
x q[4];
tdg q[4];
cx q[5], q[4];
x q[4];
tdg q[6];
cx q[6], q[7];
tdg q[7];
cx q[6], q[7];
sdg q[6];
sdg q[8];
cx q[8], q[9];
tdg q[10];
sdg q[12];
sdg q[14];
cx q[15], q[14];
tdg q[14];
sdg q[14];
tdg q[16];
cx q[18], q[19];
tdg q[19];
sdg q[18];
cx q[20], q[21];
z q[23];
tdg q[24];
z q[26];
tdg q[28];
cx q[27], q[28];
tdg q[27];
sdg q[29];
cx q[29], q[30];
sdg q[31];
z q[33];
cx q[34], q[35];
tdg q[35];
cx q[34], q[35];
sdg q[34];
tdg q[36];
cx q[39], q[38];
sdg q[38];
cx q[40], q[41];
sdg q[40];
sdg q[42];
cx q[43], q[42];
tdg q[42];
t q[61];
s q[64];
t q[82];
s q[49];
h q[61];
cx q[85], q[84];
t q[82];
cx q[84], q[63];
cx q[68], q[75];
ccx q[65], q[87], q[62];
ccx q[53], q[55], q[61];
cx q[59], q[65];
t q[85];
ccx q[59], q[68], q[65];
s q[46];
cx q[82], q[87];
h q[59];
s q[89];
t q[58];
t q[88];
cx q[63], q[66];
h q[63];
cx q[77], q[74];
cx q[64], q[48];
h q[55];
cx q[63], q[73];
t q[82];
t q[85];
ccx q[65], q[51], q[66];
ccx q[54], q[66], q[57];
t q[56];
t q[67];
ccx q[78], q[59], q[71];
cx q[76], q[72];
h q[60];
ccx q[49], q[88], q[47];
ccx q[89], q[50], q[72];
cx q[86], q[83];
ccx q[59], q[82], q[84];
t q[69];
t q[63];
s q[87];
t q[46];
t q[83];
ccx q[68], q[75], q[67];
cx q[90], q[10];
cx q[91], q[19];
cx q[92], q[51];
cx q[93], q[32];
cx q[94], q[4];
cx q[95], q[83];
cx q[96], q[41];
cx q[97], q[30];
cx q[98], q[1];
cx q[99], q[44];
cx q[100], q[6];
cx q[101], q[35];
cx q[102], q[14];
cx q[103], q[86];
cx q[104], q[69];
cx q[105], q[57];
