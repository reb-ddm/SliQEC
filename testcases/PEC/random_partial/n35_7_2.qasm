OPENQASM 2.0;
include "qelib1.inc";
qreg q[38];
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
t q[20];
t q[31];
ccx q[30], q[1], q[6];
cx q[23], q[15];
cx q[31], q[8];
cx q[28], q[23];
h q[3];
s q[10];
s q[34];
cx q[6], q[30];
s q[32];
h q[34];
ccx q[10], q[3], q[2];
cx q[23], q[10];
s q[6];
h q[8];
t q[16];
cx q[3], q[20];
t q[18];
s q[33];
ccx q[14], q[22], q[25];
t q[30];
cx q[13], q[17];
s q[34];
s q[31];
h q[5];
t q[22];
s q[18];
t q[27];
ccx q[31], q[7], q[26];
t q[27];
cx q[6], q[15];
t q[16];
s q[29];
s q[31];
ccx q[4], q[11], q[8];
h q[30];
cx q[11], q[31];
h q[11];
cx q[31], q[23];
h q[13];
cx q[21], q[19];
s q[8];
cx q[19], q[13];
cx q[27], q[9];
ccx q[8], q[33], q[4];
h q[33];
h q[10];
s q[29];
h q[3];
h q[12];
ccx q[17], q[5], q[22];
s q[17];
ccx q[1], q[18], q[13];
t q[32];
ccx q[32], q[20], q[18];
cx q[19], q[31];
cx q[29], q[31];
s q[11];
t q[3];
h q[28];
s q[26];
cx q[13], q[9];
t q[10];
s q[12];
ccx q[13], q[33], q[22];
s q[20];
ccx q[31], q[21], q[13];
cx q[27], q[31];
t q[4];
ccx q[10], q[34], q[11];
s q[9];
ccx q[33], q[28], q[6];
t q[24];
cx q[12], q[25];
s q[33];
cx q[25], q[7];
ccx q[14], q[15], q[21];
h q[14];
h q[21];
ccx q[21], q[26], q[30];
t q[18];
s q[33];
ccx q[14], q[13], q[30];
h q[23];
ccx q[7], q[16], q[22];
cx q[28], q[6];
ccx q[19], q[25], q[34];
ccx q[14], q[16], q[11];
h q[25];
h q[23];
cx q[8], q[32];
cx q[33], q[21];
s q[21];
t q[6];
h q[14];
s q[30];
t q[5];
h q[24];
cx q[33], q[2];
h q[26];
ccx q[26], q[9], q[10];
ccx q[27], q[30], q[7];
ccx q[22], q[4], q[27];
ccx q[34], q[15], q[2];
sdg q[0];
cx q[1], q[0];
tdg q[0];
cx q[3], q[2];
tdg q[2];
cx q[3], q[2];
tdg q[2];
sdg q[4];
cx q[4], q[5];
cx q[6], q[7];
tdg q[7];
tdg q[6];
cx q[6], q[7];
sdg q[8];
cx q[8], q[9];
sdg q[9];
cx q[8], q[9];
tdg q[8];
y q[10];
z q[11];
z q[13];
tdg q[14];
t q[24];
t q[30];
t q[28];
h q[24];
t q[22];
ccx q[25], q[28], q[17];
t q[31];
t q[26];
t q[27];
ccx q[32], q[21], q[34];
t q[33];
s q[31];
cx q[31], q[17];
s q[25];
t q[33];
t q[26];
t q[30];
h q[26];
cx q[35], q[24];
cx q[36], q[19];
cx q[37], q[27];
