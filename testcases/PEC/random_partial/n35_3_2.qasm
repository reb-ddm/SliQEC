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
cx q[14], q[15];
cx q[25], q[16];
s q[21];
cx q[19], q[8];
s q[8];
s q[22];
t q[9];
h q[32];
ccx q[31], q[13], q[8];
t q[9];
ccx q[22], q[32], q[17];
ccx q[33], q[10], q[6];
s q[13];
h q[27];
h q[8];
h q[3];
t q[6];
t q[14];
ccx q[26], q[21], q[6];
ccx q[8], q[15], q[6];
t q[15];
ccx q[9], q[15], q[14];
h q[8];
ccx q[18], q[8], q[23];
cx q[5], q[2];
h q[19];
t q[22];
h q[1];
s q[22];
ccx q[32], q[23], q[0];
s q[12];
s q[24];
t q[9];
ccx q[30], q[27], q[15];
h q[15];
s q[27];
t q[10];
h q[17];
s q[17];
h q[24];
t q[7];
cx q[18], q[19];
s q[32];
ccx q[11], q[34], q[0];
h q[17];
h q[14];
s q[19];
s q[17];
s q[11];
t q[31];
cx q[14], q[2];
t q[28];
ccx q[23], q[24], q[12];
cx q[20], q[14];
ccx q[20], q[29], q[8];
h q[14];
h q[13];
ccx q[16], q[5], q[23];
t q[15];
t q[30];
s q[6];
s q[3];
h q[5];
t q[17];
cx q[3], q[5];
h q[8];
h q[18];
ccx q[11], q[16], q[1];
h q[19];
ccx q[0], q[1], q[21];
s q[1];
ccx q[8], q[21], q[9];
s q[1];
h q[19];
cx q[24], q[27];
t q[11];
h q[0];
h q[14];
cx q[30], q[29];
s q[27];
h q[6];
ccx q[34], q[22], q[29];
h q[17];
cx q[19], q[4];
ccx q[9], q[19], q[21];
t q[18];
cx q[3], q[9];
cx q[9], q[16];
cx q[11], q[23];
h q[4];
h q[22];
ccx q[34], q[11], q[0];
h q[23];
ccx q[27], q[31], q[9];
t q[33];
cx q[16], q[18];
h q[15];
h q[8];
t q[14];
s q[23];
t q[7];
t q[7];
ccx q[29], q[21], q[32];
t q[24];
ccx q[8], q[10], q[15];
tdg q[1];
sdg q[0];
cx q[0], q[1];
tdg q[2];
cx q[3], q[2];
sdg q[3];
tdg q[2];
sdg q[5];
sdg q[4];
cx q[5], q[4];
sdg q[4];
sdg q[6];
tdg q[6];
cx q[9], q[8];
tdg q[8];
sdg q[9];
cx q[9], q[8];
tdg q[8];
sdg q[10];
cx q[12], q[13];
sdg q[12];
x q[14];
cx q[16], q[15];
tdg q[15];
cx q[16], q[15];
sdg q[16];
sdg q[15];
h q[21];
h q[19];
t q[17];
cx q[28], q[17];
ccx q[24], q[22], q[33];
ccx q[22], q[25], q[30];
cx q[18], q[26];
h q[21];
t q[28];
ccx q[29], q[32], q[31];
ccx q[21], q[29], q[20];
t q[32];
ccx q[32], q[34], q[19];
t q[24];
cx q[22], q[29];
cx q[20], q[28];
cx q[17], q[20];
ccx q[29], q[27], q[24];
cx q[35], q[17];
cx q[36], q[31];
cx q[37], q[10];
cx q[38], q[21];
cx q[39], q[0];
cx q[40], q[9];
