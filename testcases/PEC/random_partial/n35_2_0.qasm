OPENQASM 2.0;
include "qelib1.inc";
qreg q[35];
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
s q[25];
s q[21];
s q[32];
cx q[24], q[9];
ccx q[0], q[29], q[22];
ccx q[8], q[25], q[4];
s q[5];
t q[16];
t q[32];
t q[1];
ccx q[28], q[17], q[31];
s q[30];
h q[32];
s q[34];
cx q[7], q[20];
ccx q[7], q[16], q[21];
h q[0];
cx q[22], q[21];
s q[29];
ccx q[26], q[10], q[4];
h q[23];
h q[10];
h q[16];
s q[14];
ccx q[33], q[27], q[2];
t q[32];
ccx q[30], q[19], q[8];
t q[12];
cx q[12], q[30];
ccx q[22], q[0], q[18];
t q[16];
ccx q[29], q[24], q[6];
cx q[6], q[7];
s q[24];
s q[30];
t q[28];
t q[18];
h q[26];
t q[14];
cx q[18], q[28];
ccx q[13], q[14], q[33];
ccx q[6], q[17], q[8];
h q[21];
h q[32];
t q[31];
t q[7];
t q[4];
cx q[7], q[19];
cx q[29], q[18];
s q[23];
ccx q[5], q[22], q[6];
cx q[29], q[26];
s q[22];
s q[11];
cx q[13], q[9];
h q[4];
ccx q[11], q[22], q[0];
ccx q[31], q[7], q[2];
h q[0];
t q[4];
ccx q[8], q[14], q[26];
cx q[0], q[1];
h q[31];
s q[22];
s q[7];
ccx q[24], q[8], q[26];
ccx q[5], q[10], q[1];
s q[9];
s q[33];
t q[9];
t q[23];
s q[7];
ccx q[7], q[28], q[19];
s q[16];
cx q[8], q[30];
cx q[1], q[7];
h q[29];
cx q[30], q[7];
s q[21];
ccx q[29], q[31], q[9];
h q[10];
s q[5];
h q[0];
t q[19];
ccx q[8], q[5], q[12];
t q[19];
h q[7];
t q[27];
ccx q[17], q[29], q[7];
s q[10];
ccx q[33], q[14], q[2];
ccx q[16], q[24], q[15];
s q[30];
cx q[32], q[27];
s q[30];
s q[22];
h q[15];
t q[15];
ccx q[21], q[2], q[26];
h q[28];
s q[13];
h q[12];
cx q[33], q[10];
s q[30];
t q[13];
cx q[0], q[1];
sdg q[1];
tdg q[1];
cx q[0], q[1];
sdg q[0];
sdg q[2];
cx q[3], q[2];
tdg q[2];
cx q[3], q[2];
tdg q[2];
tdg q[5];
tdg q[4];
cx q[4], q[5];
cx q[7], q[6];
tdg q[6];
sdg q[6];
sdg q[9];
cx q[8], q[9];
sdg q[9];
cx q[8], q[9];
tdg q[8];
cx q[10], q[11];
sdg q[10];
x q[12];
cx q[14], q[13];
sdg q[13];
cx q[14], q[13];
tdg q[13];
tdg q[16];
cx q[15], q[16];
tdg q[16];
tdg q[15];
t q[23];
h q[17];
t q[34];
t q[27];
s q[19];
t q[27];
h q[26];
s q[34];
cx q[31], q[22];
cx q[22], q[31];
cx q[30], q[27];
ccx q[34], q[32], q[17];
h q[21];
h q[26];
cx q[17], q[32];
ccx q[29], q[32], q[22];
h q[24];
ccx q[30], q[20], q[34];
