OPENQASM 2.0;
include "qelib1.inc";
qreg q[59];
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
s q[2];
cx q[34], q[43];
cx q[35], q[40];
h q[48];
s q[48];
s q[25];
t q[39];
ccx q[5], q[36], q[1];
s q[1];
ccx q[15], q[26], q[7];
ccx q[12], q[48], q[8];
cx q[30], q[2];
s q[34];
ccx q[9], q[25], q[7];
ccx q[29], q[30], q[21];
cx q[2], q[24];
h q[36];
s q[16];
ccx q[3], q[22], q[30];
h q[34];
s q[34];
cx q[35], q[22];
cx q[27], q[7];
h q[41];
s q[29];
s q[15];
t q[9];
h q[49];
h q[40];
ccx q[23], q[48], q[35];
s q[41];
cx q[24], q[44];
h q[49];
h q[39];
cx q[29], q[49];
ccx q[47], q[12], q[39];
s q[27];
h q[43];
ccx q[44], q[33], q[22];
h q[48];
ccx q[8], q[39], q[9];
cx q[37], q[36];
h q[48];
h q[30];
t q[7];
ccx q[20], q[17], q[22];
s q[19];
s q[38];
s q[24];
cx q[7], q[5];
cx q[48], q[2];
cx q[7], q[27];
h q[32];
s q[11];
ccx q[43], q[7], q[2];
ccx q[18], q[27], q[29];
s q[17];
t q[8];
h q[12];
h q[1];
ccx q[38], q[46], q[12];
t q[22];
cx q[29], q[10];
cx q[30], q[46];
s q[26];
ccx q[8], q[37], q[30];
h q[36];
t q[3];
ccx q[16], q[23], q[45];
h q[15];
ccx q[1], q[14], q[22];
ccx q[28], q[13], q[38];
cx q[36], q[26];
ccx q[3], q[21], q[44];
t q[48];
t q[38];
cx q[45], q[29];
s q[8];
t q[11];
cx q[9], q[29];
s q[47];
h q[27];
cx q[42], q[0];
t q[38];
t q[34];
t q[9];
ccx q[6], q[18], q[15];
s q[23];
ccx q[45], q[34], q[39];
ccx q[38], q[27], q[12];
h q[10];
ccx q[11], q[10], q[14];
cx q[47], q[22];
ccx q[49], q[1], q[22];
s q[13];
s q[37];
s q[12];
cx q[43], q[9];
s q[37];
h q[33];
s q[3];
t q[18];
cx q[40], q[49];
cx q[35], q[49];
t q[12];
h q[46];
t q[12];
ccx q[17], q[27], q[1];
cx q[29], q[25];
s q[45];
h q[29];
ccx q[45], q[40], q[11];
ccx q[31], q[8], q[45];
ccx q[27], q[26], q[22];
cx q[31], q[44];
ccx q[16], q[48], q[13];
s q[16];
t q[7];
t q[4];
h q[49];
ccx q[26], q[15], q[24];
t q[35];
s q[15];
s q[31];
cx q[49], q[39];
t q[44];
cx q[19], q[16];
cx q[37], q[33];
cx q[23], q[31];
ccx q[30], q[1], q[16];
ccx q[2], q[17], q[32];
s q[14];
ccx q[15], q[2], q[4];
ccx q[5], q[23], q[1];
s q[36];
s q[24];
ccx q[21], q[19], q[43];
s q[30];
h q[22];
h q[27];
ccx q[8], q[18], q[11];
h q[6];
s q[19];
cx q[47], q[18];
cx q[22], q[4];
t q[0];
cx q[15], q[39];
t q[17];
cx q[3], q[41];
h q[48];
y q[0];
cx q[2], q[3];
tdg q[3];
sdg q[2];
tdg q[2];
cx q[2], q[3];
tdg q[4];
cx q[7], q[8];
sdg q[7];
cx q[10], q[9];
tdg q[10];
tdg q[9];
cx q[10], q[9];
tdg q[9];
cx q[12], q[13];
sdg q[12];
y q[15];
cx q[16], q[17];
sdg q[16];
sdg q[17];
cx q[16], q[17];
tdg q[16];
tdg q[19];
cx q[18], q[19];
tdg q[19];
cx q[18], q[19];
sdg q[18];
cx q[20], q[21];
x q[22];
h q[26];
h q[25];
cx q[44], q[30];
h q[48];
cx q[30], q[45];
cx q[40], q[30];
t q[29];
s q[49];
h q[49];
t q[30];
s q[27];
cx q[41], q[47];
ccx q[26], q[48], q[40];
s q[26];
cx q[35], q[28];
cx q[28], q[26];
t q[26];
t q[36];
cx q[30], q[33];
t q[36];
h q[32];
cx q[39], q[25];
s q[27];
t q[45];
ccx q[34], q[43], q[32];
cx q[50], q[2];
cx q[51], q[21];
cx q[52], q[30];
cx q[53], q[7];
cx q[54], q[3];
cx q[55], q[19];
cx q[56], q[38];
cx q[57], q[44];
cx q[58], q[36];
