## Heat distribution test
This test is an example of a stencil-based problem. It implements
a numerical algorithm that solves PDEs simulating heat distribution process.

## How to run

1. sgc heat.sol
2. mpicxx out.cpp <options...> -o test
3. g++ png\_writer.cpp -lpng -p png\_writer
4. ./png\_writer 10 10 100 100

The default command-line arguments for png converter:

`./png <nb> <mb> <blocn_n> <block_m>`, where `nb`, `mb`, `block_n`, `block_m`
are the names of the corresponding variables declared in heat.sol.
