OPENQASM 2.0;
include "qelib1.inc";
qreg q[41];
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
cx q[20], q[8];
cx q[15], q[31];
h q[15];
h q[21];
cx q[17], q[20];
h q[7];
t q[15];
s q[4];
ccx q[1], q[18], q[14];
t q[1];
h q[25];
t q[4];
cx q[10], q[22];
ccx q[30], q[16], q[21];
h q[28];
s q[18];
s q[30];
cx q[17], q[12];
ccx q[7], q[17], q[0];
ccx q[11], q[9], q[0];
ccx q[13], q[3], q[15];
s q[19];
s q[16];
ccx q[27], q[0], q[7];
h q[19];
t q[22];
h q[33];
h q[16];
cx q[30], q[10];
cx q[2], q[23];
cx q[11], q[28];
ccx q[23], q[34], q[17];
cx q[6], q[9];
s q[29];
s q[9];
h q[0];
t q[22];
h q[1];
cx q[27], q[19];
t q[17];
cx q[19], q[7];
ccx q[30], q[18], q[21];
t q[11];
h q[1];
cx q[28], q[1];
h q[10];
cx q[28], q[8];
cx q[25], q[26];
cx q[9], q[4];
t q[8];
t q[20];
h q[33];
t q[13];
t q[29];
s q[26];
t q[19];
cx q[24], q[27];
s q[32];
cx q[30], q[6];
h q[24];
t q[1];
cx q[34], q[28];
h q[32];
h q[13];
ccx q[13], q[18], q[15];
cx q[26], q[8];
cx q[14], q[33];
cx q[2], q[5];
h q[23];
s q[29];
h q[4];
t q[17];
s q[21];
ccx q[23], q[24], q[4];
cx q[18], q[25];
ccx q[25], q[2], q[0];
ccx q[24], q[9], q[13];
ccx q[19], q[7], q[15];
t q[3];
ccx q[25], q[1], q[21];
t q[17];
cx q[16], q[34];
ccx q[17], q[23], q[2];
t q[5];
s q[8];
cx q[16], q[30];
t q[21];
ccx q[25], q[32], q[16];
h q[15];
h q[15];
h q[19];
s q[19];
cx q[1], q[17];
t q[12];
ccx q[7], q[16], q[8];
t q[9];
cx q[21], q[28];
t q[18];
h q[16];
s q[17];
h q[32];
ccx q[15], q[2], q[19];
cx q[23], q[12];
s q[8];
cx q[9], q[25];
cx q[0], q[1];
sdg q[1];
cx q[0], q[1];
sdg q[0];
cx q[3], q[2];
sdg q[2];
sdg q[4];
tdg q[4];
x q[7];
tdg q[9];
cx q[8], q[9];
sdg q[9];
sdg q[8];
cx q[11], q[10];
tdg q[10];
cx q[11], q[10];
sdg q[10];
tdg q[10];
tdg q[13];
cx q[12], q[13];
sdg q[13];
cx q[12], q[13];
tdg q[12];
sdg q[14];
tdg q[14];
cx q[15], q[14];
sdg q[14];
s q[28];
s q[33];
h q[27];
cx q[17], q[21];
cx q[29], q[32];
cx q[28], q[34];
t q[21];
t q[34];
ccx q[18], q[17], q[31];
ccx q[33], q[30], q[27];
cx q[20], q[32];
ccx q[26], q[29], q[18];
s q[22];
cx q[31], q[21];
h q[30];
t q[21];
ccx q[25], q[24], q[21];
s q[28];
cx q[35], q[20];
cx q[36], q[2];
cx q[37], q[8];
cx q[38], q[10];
cx q[39], q[17];
cx q[40], q[21];
