// `commands.cpp` - An attempt to make a C++ version of woflang

#include <iostream>
#include <unordered_map>
#include <functional>
#include <cmath>
#include <bitset>
#include <cstring>
#include <algorithm>
#include "utils.h"
#include "commands.h"

// Function prototypes
void cmdCmp(double* operands, int numOperands);
void cmdAdd(double* operands, int numOperands);
void cmdSub(double* operands, int numOperands);
void cmdMul(double* operands, int numOperands);
void cmdDiv(double* operands, int numOperands);
void cmdEq(double* operands, int numOperands);
void cmdDNE(double* operands, int numOperands);
void cmdGT(double* operands, int numOperands);
void cmdLT(double* operands, int numOperands);
void cmdGTE(double* operands, int numOperands);
void cmdLTE(double* operands, int numOperands);
void cmdAnd(double* operands, int numOperands);
void cmdOr(double* operands, int numOperands);
void cmdXor(double* operands, int numOperands);
void cmdNot(double* operands, int numOperands);
void cmdNand(double* operands, int numOperands);
void cmdNor(double* operands, int numOperands);
void cmdIf(double* operands, int numOperands);
void cmdThen(double* operands, int numOperands);
void cmdElse(double* operands, int numOperands);
void cmdWhile(double* operands, int numOperands);
void cmdFor(int start, int end, int step, std::function<void(int)> body);
void cmdBuf(double* operands, int numOperands);
void cmdInhibit(double* operands, int numOperands);
void cmdEquiv(double* operands, int numOperands);
void cmdId(double* operands, int numOperands);
void cmdBin2Fuz(double* operands, int numOperands);
void cmdFuz2Bin(double* operands, int numOperands);
void cmdBool(double* operands, int numOperands);
void cmdCount(double* operands, int numOperands);
void cmdStr(double* operands, int numOperands);
void cmdVectAdd(double* operands, int numOperands);
void cmdVectSub(double* operands, int numOperands);
void cmdDotProd(double* operands, int numOperands);
void cmdMtrxDet(double* operands, int numOperands);
void cmdMtrxInv(double* operands, int numOperands);
void cmd3x3Det(double* operands, int numOperands);
void cmdMtrxAnd(double* operands, int numOperands);
void cmdBreak(double* operands, int numOperands);
void cmdEnd(double* operands, int numOperands);
void cmdExit(double* operands, int numOperands);
void cmdPrompt(double* operands, int numOperands);
void cmdDone(double* operands, int numOperands);
void cmdProxy(double* operands, int numOperands);
void cmdGetProxy(double* operands, int numOperands);
void cmdServer(double* operands, int numOperands);
void cmdUserAgent(double* operands, int numOperands);
void cmdFnct(double* operands, int numOperands);
void cmdAttribute(double* operands, int numOperands);
void cmdType(double* operands, int numOperands);
void cmdWait(double* operands, int numOperands);
void cmdSleep(double* operands, int numOperands);
void cmdPrint(double* operands, int numOperands);
void cmdChar(double* operands, int numOperands);
void cmdRead(double* operands, int numOperands);
void cmdMove(double* operands, int numOperands);
void cmdNull(double* operands, int numOperands);
void cmdPath(double* operands, int numOperands);

// Define any conditionally compiled code
#ifdef ENABLE_CORE_MATH
void cmdCoreMath(double* operands, int numOperands) {
    if (operands == nullptr || numOperands < 1) {
        reportError(ERR_MEM_INVALID_OP, "cmdCoreMath - Invalid operands or number of operands");
        return;
    }

// Utility function to clamp values to the range [0, 100]
double clamp(double value) {
    return std::max(0.0, std::min(100.0, value));
}

// Command functions
void cmdCmp(double* operands, int numOperands) {
    if (numOperands >= 2) {
        double result = (operands[0] > operands[1]) ? 100 : 0;  // Fuzzy comparison
        std::cout << "Result (CMP): " << result << std::endl;
    }
}

void cmdAdd(double* operands, int numOperands) {
    if (numOperands >= 2) {
        double result = operands[0] + operands[1];
        result = std::min(result, 100.0);  // Clamp to max 100
        std::cout << "Result (ADD): " << result << std::endl;
    }
}

void cmdSub(double* operands, int numOperands) {
    if (numOperands >= 2) {
        double result = operands[0] - operands[1];
        result = std::max(result, 0.0);  // Clamp to min 0
        std::cout << "Result (SUBTRACT): " << result << std::endl;
    }
}

void cmdMul(double* operands, int numOperands) {
    if (numOperands >= 2) {
        double result = (operands[0] * operands[1]) / 100.0;  // Normalize product
        std::cout << "Result (MULTIPLY): " << result << std::endl;
    }
}

#include "err_chk.h" // Include the error checking header

void cmdDiv(double* operands, int numOperands) {
    if (numOperands >= 2) {
        if (operands[1] == 0) {
            reportError(ERR_MATH_DIV_ZERO, "cmdDiv"); // Report division by zero error
            return;
        }
        double result = (operands[0] / operands[1]) * 100;  // Normalize division
        std::cout << "Result (DIVIDE): " << result << std::endl;
    } else {
        reportError(ERR_UNRECOGNIZED, "cmdDiv - Not enough operands");
    }
}

void cmdEq(double* operands, int numOperands) {
    if (numOperands >= 2) {
        double result = (operands[0] == operands[1]) ? 100 : 0;  // 100 if equal
        std::cout << "Result (EQUALS): " << result << std::endl;
    }
}

void cmdDNE(double* operands, int numOperands) {
    if (numOperands >= 2) {
        double result = (operands[0] != operands[1]) ? 100 : 0;  // 100 if not equal
        std::cout << "Result (DOES NOT EQUAL): " << result << std::endl;
    }
}

void cmdGT(double* operands, int numOperands) {
    if (numOperands >= 2) {
        double result = (operands[0] > operands[1]) ? 100 : 0;  // 100 if greater
        std::cout << "Result (GREATER THAN): " << result << std::endl;
    }
}

void cmdLT(double* operands, int numOperands) {
    if (numOperands >= 2) {
        double result = (operands[0] < operands[1]) ? 100 : 0;  // 100 if less
        std::cout << "Result (LESS THAN): " << result << std::endl;
    }
}

void cmdGTE(double* operands, int numOperands) {
    if (numOperands >= 2) {
        double result = (operands[0] >= operands[1]) ? 100 : 0;  // 100 if greater or equal
        std::cout << "Result (GREATER THAN OR EQUAL): " << result << std::endl;
    }
}

void cmdLTE(double* operands, int numOperands) {
    if (numOperands >= 2) {
        double result = (operands[0] <= operands[1]) ? 100 : 0;  // 100 if less or equal
        std::cout << "Result (LESS THAN OR EQUAL): " << result << std::endl;
    }
}

void cmdAnd(double* operands, int numOperands) {
    if (numOperands >= 2) {
        double result = std::min(operands[0], operands[1]);  // Fuzzy AND: MIN function
        std::cout << "Result (AND): " << result << std::endl;
    }
}

void cmdOr(double* operands, int numOperands) {
    if (numOperands >= 2) {
        double result = std::max(operands[0], operands[1]);  // Fuzzy OR: MAX function
        std::cout << "Result (OR): " << result << std::endl;
    }
}

void cmdXor(double* operands, int numOperands) {
    if (numOperands >= 2) {
        double result = std::fabs(operands[0] - operands[1]);  // Fuzzy XOR: Absolute difference
        std::cout << "Result (XOR): " << result << std::endl;
    }
}

void cmdNot(double* operands, int numOperands) {
    if (numOperands >= 1) {
        double result = 100 - operands[0];  // Fuzzy NOT: Inversion
        std::cout << "Result (NOT): " << result << std::endl;
    }
}

void cmdNand(double* operands, int numOperands) {
    if (numOperands >= 2) {
        double andResult = std::min(operands[0], operands[1]);
        double result = 100 - andResult;  // Fuzzy NAND: NOT of AND
        std::cout << "Result (NAND): " << result << std::endl;
    }
}

void cmdNor(double* operands, int numOperands) {
    if (numOperands >= 2) {
        double orResult = std::max(operands[0], operands[1]);
        double result = 100 - orResult;  // Fuzzy NOR: NOT of OR
        std::cout << "Result (NOR): " << result << std::endl;
    }
}

void cmdIf(double* operands, int numOperands) {
    if (numOperands >= 3) {
        double result = (operands[0] > 0) ? operands[1] : operands[2];  // Fuzzy IF: Conditional selection
        std::cout << "Result (IF): " << result << std::endl;
    }
}

void cmdThen(double* operands, int numOperands) {
    if (numOperands >= 2) {
        std::cout << "Result (THEN): " << clamp(operands[0] ? operands[1] : 0) << std::endl;
    }
}

void cmdElse(double* operands, int numOperands) {
    if (numOperands >= 1) {
        std::cout << "Result (ELSE): " << clamp(operands[0]) << std::endl;
    }
}

void cmdWhile(double* operands, int numOperands) {
    if (numOperands >= 2) {
        while (operands[0] > 0) {
            std::cout << "Executing (WHILE) with operand: " << clamp(operands[0]) << std::endl;
            operands[0]--;  // Example loop decrement
        }
    }
}

void cmdFor(int start, int end, int step, std::function<void(int)> body) {
    for (int i = start; i < end; i += step) {
        body(i);
    }
}

void cmdBuf(double* operands, int numOperands) {
    if (numOperands >= 2) {
        std::cout << "Result (BUFFER): " << clamp((operands[0] = operands[1])) << std::endl;
    }
}

void cmdInhibit(double* operands, int numOperands) {
    if (numOperands >= 2) {
        std::cout << "Result (INHIBIT): " << clamp((operands[0] = operands[0])) << std::endl;
    }
}

void cmdEquiv(double* operands, int numOperands) {
    if (numOperands >= 2) {
        std::cout << "Result (EQUIV): " << clamp(operands[0] == operands[1] ? 100 : 0) << std::endl;
    }
}

void cmdId(double* operands, int numOperands) {
    if (numOperands >= 1) {
        std::cout << "Result (IDENTITY): " << clamp(operands[0]) << std::endl;
    }
}

void cmdBool(double* operands, int numOperands) {
    if (numOperands >= 1) {
        std::cout << "Result (BOOLEAN): " << (operands[0] != 0 ? "true" : "false") << std::endl;
    }
}

void cmdBin2Fuz(double* operands, int numOperands) {
    if (numOperands >= 1) {
        std::cout << "Result (BINARY-2-FUZZY): " << clamp(100 - operands[0]) << std::endl;
    }
}

void cmdFuz2Bin(double* operands, int numOperands) {
    if (numOperands >= 1) {
        std::cout << "Result (FUZZY-2-BINARY): " << clamp(operands[0]) << std::endl;  // Assuming fuzzy value is already in [0, 100]
    }
}

void cmdVectAdd(double* operands, int numOperands) {
    if (numOperands >= 2) {
        std::cout << "Result (VECTOR ADD): " << clamp(operands[0] + operands[1]) << std::endl;
    }
}

void cmdVectSub(double* operands, int numOperands) {
    if (numOperands >= 2) {
        std::cout << "Result (VECTOR SUB): " << clamp(operands[0] - operands[1]) << std::endl;
    }
}

void cmdDotProd(double* operands, int numOperands) {
    if (numOperands >= 2) {
        std::cout << "Result (DOT PRODUCT): " << clamp(operands[0] * operands[1]) << std::endl;
    }
}

void cmdMtrxDet(double* operands, int numOperands) {
    if (numOperands >= 2) {
        // Example operation, replace with actual matrix determinant calculation
        std::cout << "Result (MATRIX DET): " << clamp(operands[0] * operands[1]) << std::endl;
    }
}

void cmdMtrxInv(double* operands, int numOperands) {
    if (numOperands >= 2) {
        // Example operation, replace with actual matrix inversion calculation
        std::cout << "Result (MATRIX INV): " << clamp(operands[0] / operands[1]) << std::endl;
    }
}

void cmd3x3Det(double* operands, int numOperands) {
    if (numOperands >= 2) {
        // Example operation, replace with actual 3x3 matrix determinant calculation
        std::cout << "Result (3X3 DET): " << clamp(operands[0] * operands[1]) << std::endl;
    }
}

void cmdMtrxAnd(double* operands, int numOperands) {
    if (numOperands >= 2) {
        // Example operation, replace with actual matrix AND calculation
        std::cout << "Result (MATRIX AND): " << clamp(operands[0] && operands[1] ? 100 : 0) << std::endl;
    }
}

void cmdBreak(double* operands, int numOperands) {
    if (numOperands >= 2) {
        std::cout << "Result (BREAK): " << clamp(operands[0]) << std::endl;
    }
}

void cmdEnd(double* operands, int numOperands) {
    if (numOperands >= 2) {
        std::cout << "Result (END): " << clamp(operands[0]) << std::endl;
    }
}

void cmdExit(double* operands, int numOperands) {
    if (numOperands >= 2) {
        std::cout << "Result (EXIT): " << clamp(operands[0]) << std::endl;
    }
}

void cmdPrompt(double* operands, int numOperands) {
    if (numOperands >= 2) {
        std::cout << "Result (PROMPT): " << clamp(operands[0]) << std::endl;
    }
}

void cmdDone(double* operands, int numOperands) {
    if (numOperands >= 2) {
        std::cout << "Result (DONE): " << clamp(operands[0]) << std::endl;
    }
}


void cmdProxy(double* operands, int numOperands) {
    if (numOperands >= 2) {
        std::cout << "Result (PROXY): " << clamp(operands[0]) << std::endl;
    }
}

void cmdGetProxy(double* operands, int numOperands) {
    if (numOperands >= 2) {
        std::cout << "Result (GET PROXY): " << clamp(operands[0]) << std::endl;
    }
}

void cmdServer(double* operands, int numOperands) {
    if (numOperands >= 2) {
        std::cout << "Result (SERVER): " << clamp(operands[0]) << std::endl;
    }
}

void cmdUserAgent(double* operands, int numOperands) {
    if (numOperands >= 2) {
        std::cout << "Result (USER AGENT): " << clamp(operands[0]) << std::endl;
    }
}

void cmdFnct(double* operands, int numOperands) {
    if (numOperands >= 2) {
        std::cout << "Result (FUNCTION): " << clamp(operands[0]) << std::endl;
    }
}

void cmdAttribute(double* operands, int numOperands) {
    if (numOperands >= 2) {
        std::cout << "Result (ATTRIBUTE): " << clamp(operands[0]) << std::endl;
    }
}

void cmdType(double* operands, int numOperands) {
    if (numOperands >= 2) {
        std::cout << "Result (TYPE): " << clamp(operands[0]) << std::endl;
    }
}

void cmdWait(double* operands, int numOperands) {
    if (numOperands >= 2) {
        std::cout << "Result (WAIT): " << clamp(operands[0]) << std::endl;
    }
}

void cmdSleep(double* operands, int numOperands) {
    if (numOperands >= 2) {
        std::cout << "Result (SLEEP): " << clamp(operands[0]) << std::endl;
    }
}

void cmdPrint(double* operands, int numOperands) {
    if (numOperands >= 2) {
        std::cout << "Result (PRINT): " << clamp(operands[0]) << std::endl;
    }
}

void cmdChar(double* operands, int numOperands) {
    if (numOperands >= 2) {
        std::cout << "Result (CHAR): " << clamp(operands[0]) << std::endl;
    }
}

void cmdRead(double* operands, int numOperands) {
    if (numOperands >= 2) {
        std::cout << "Result (READ): " << clamp(operands[0]) << std::endl;
    }
}

void cmdMove(double* operands, int numOperands) {
    if (numOperands >= 2) {
        std::cout << "Result (MOVE): " << clamp(operands[0]) << std::endl;
    }
}

void cmdNull(double* operands, int numOperands) {
    if (numOperands >= 2) {
        std::cout << "Result (NULL): " << clamp(operands[0]) << std::endl;
    }
}

void cmdPath(double* operands, int numOperands) {
    if (numOperands >= 2) {
        std::cout << "Result (PATH): " << clamp(operands[0]) << std::endl;
    }
}

void executeCommand(CommandType commandType, double* operands, int numOperands) {
    auto it = commandMap.find(commandType);
    if (it != commandMap.end()) 
        it->second(operands, numOperands);
    else {
        std::cerr << "Unknown command" << std::endl;
        std::cerr << "Error code: " << ERR_UNRECOGNIZED << std::endl;
    }
}

int main() {
    double operands[] = {5.0, 3.0};
    int numOperands = 2;

    executeCommand(0x01, operands, numOperands); // Example usage for ADD
    executeCommand(0x04, operands, numOperands);

    return 0;
    } // Ensure this closing brace is present
}
#endif

#ifdef ENABLE_CORE_MATH
void cmdCoreMath(double* operands, int numOperands) {
// Core math operations below:
}
std::unordered_map<CommandType, CommandFunction> commandMap = {
    {CMD_CMP, cmdCmp},
    {0x2F60, cmdCmp},  // ⽠ (melon) symbol for cmp
    {CMD_ADD, cmdAdd},
    {0x2FD1, cmdAdd},  // 齊 (even) symbol for add
    {CMD_SUB, cmdSub},
    {0x2FAA, cmdSub},  // 隶 (slave) symbol for subtract
    {CMD_MUL, cmdMul},
    {0x2F02, cmdMul},  // 丶 (dot) symbol for multiply
    {CMD_DIV, cmdDiv},
    {0x2F03, cmdDiv},  // 丿 (slash) symbol for divide
    {CMD_GT, cmdGT},       
    {0x2F24, cmdGT},  // 大 (big) symbol for greater than
    {CMD_LT, cmdLT},          
    {0x2F29, cmdLT},  // 小 (small) symbol for less than
    {CMD_EQ, cmdEq},          
    {0x2F22, cmdEq},  // 卜 (divination) symbol for equals
    {CMD_DNE, cmdDNE},         
    {0x2F18, cmdDNE},  // 丨 (line) symbol for does not equal
    {CMD_GTE, cmdGTE},        
    {0x2F9B, cmdGTE},  // 走 (run) symbol for greater than or equal
    {CMD_LTE, cmdLTE},        
    {0x2F21, cmdLTE},  // 夊 (slowly) symbol for less than or equal
    {CMD_AND, cmdAnd},         
    {0x2F7D, cmdAnd},  // 而 (and) symbol for AND
    {CMD_OR, cmdOr},          
    {0x6216, cmdOr},  // 或 (or) symbol for OR
    {CMD_XOR, cmdXor},         
    {0x2295, cmdXor},  // ⊕ (direct sum) symbol for XOR
    {CMD_NOT, cmdNot},         
    {0x2F46, cmdNot},  // 无 (not) symbol for NOT
    {CMD_NAND, cmdNand},        
    {0x2F46 - 0x200D - 0x2F7D, cmdNand},  // 无‍而 (not-and) symbol for NAND
    {CMD_NOR, cmdNor},         
    {0x2F46 - 0x200D - 0x6216, cmdNor},  // 无‍或 (not-or) symbol for NOR
    {CMD_IF, cmdIf},          
    {0x82E5, cmdIf},  // 若 (if) symbol for IF
    {CMD_THEN, cmdThen},        
    {0x5219, cmdThen},  // 则 (then) symbol for THEN
    {CMD_ELSE, cmdElse},        
    {0x53E6, cmdElse},  // 另 (other) symbol for ELSE
    {CMD_WHILE, cmdWhile},       
    {0x2F7C, cmdWhile},  // 老 (slowly) symbol for while
    {CMD_FOR, cmdFor},         
    {0x2F84, cmdFor},  // 至 (old) symbol for for
    {CMD_BUF, cmdBuf},         
    {0x2F87,cmdBuf},  // 舛  (arrive) symbol for buffer
    {CMD_INHIBIT, cmdInhibit},     
    {0x2F7D, cmdInhibit},  // 自 (oppose) symbol for inhibit
    {CMD_EQUIV, cmdEquiv},       
    {0x2F3C, cmdEquiv},  // 又 (self) symbol for equivalent
    {CMD_ID, cmdId},         
    {0x5165, cmdId},  // 入 (again) symbol for identity
    {CMD_BIN2FUZ, cmdBin2Fuz},     
    {0x2F4C, cmdBin2fuz},  // 止 (in) symbol for binary to fuzzy logic
    {CMD_FUZ2BIN, cmdFuz2Bin},     
    {0x51FA, cmdFuz2bin},  // 出 (stop) symbol for fuzzy to binary logic
    {CMD_BOOL, cmdBool},       
    {0x2F52/2F95, cmdB0ol},  // 氏/ 谷 (out) symbol for boolean
    {CMD_COUNT, cmdCount},       
    {0x2F30, cmdCount},  //己  (clan/valley) symbol for count
    {CMD_STR, cmdStr},         
    {0x2FA1, cmdStr},  // 辵 (oneself) symbol for string
    {CMD_VECTADD, cmdVectAdd},     
    {0x2F50, cmdVectadd},  // 比  (walk) symbol for vector aaddition
    {CMD_VECTSUB, cmdVectSub},     
    {0x7D42, cmdVectsub},  // 終 (compare) symbol for vector subtraction
    {CMD_DOTPROD, cmdDotProd},     
    {0x97FF, cmdDotprod},  // 響 (and) symbol for dot products
    {CMD_MTRXDET, cmdMtrxDet},    
    {0x2ED4, cmdMatrixdet},  // 釆 (echo) symbol for matrix determinant
    {CMD_MTRXINV, cmdMtrxInv},     
    {0x5369, cmdMatrixinv},  // 卩 (stand) symbol for matrix inversion
    {CMD_3X3DET, cmd3x3Det},      
    {0x62BC, cmd3x3det},  // 押 (seal/go out) symbol 3x3 matix precalced
    {CMD_MTRXAND, cmdMtrxAnd},     
    {0x2F40, cmdMatrixand},  //支  (push) symbol for matrix logical AND
    {CMD_BREAK, cmdBreak},      
    {0x2F70, cmdBreak},  //⽰  (branch) symbol for break
    {CMD_END, cmdEnd},        
    {0x7720, cmdEnd},  //眠  (spirit) symbol for end
    {CMD_EXIT, cmdExit},        
    {0x2E92, cmdExit},  // ⺒  (drowsy) symbol for exit
    {CMD_PROMPT, cmdPrompt},      
    {0x2EA8, cmdPrompt},  // ⺨ (snake) symbol for prompt
    {CMD_DONE, cmdDone},                
    {0x2F6A, cmdDone},  // ⽪ (dog symbol for done
    {CMD_PROXY, cmdProxy},       
    {0x2F53, cmdProxy},  // ⽓ (skin) symbol for proxy
    {CMD_GETPROXY, cmdGetProxy},    
    {0x5F85, cmdGetproxy},  //待  (steam) symbol for get proxy
    {CMD_SERVER, cmdServer},  
    {0x5EE4, cmdServer},  // ⻤ (wait) symbol for server
    {CMD_USERAGENT, cmdUserAgent},
    {0x2E99, cmdUseragent},  //⺙  (ghost) symbol for user agent
    {CMD_FNCT, cmdFnct},
    {0x2E96, cmdFnct},  // ⺆ (rap) symbol for function
    {CMD_ATTRIBUTE, cmdAttribute},
    {0x2E98, cmdAttribute},  //⺘  (box) symbol for attribute
    {CMD_TYPE, cmdType},
    {0x2F2F, cmdType},  // 工 (hand) symbol for type
    {CMD_WAIT, cmdWait},
    {0x2EA9, cmdWait},  //⺣  (work) symbol for wait
    {CMD_SLEEP, cmdSleep},
    {0x5DE1, cmdSleep},  // 巡 (fire) symbol for sleep
    {CMD_PRINT, cmdPrint},
    {0x2FAE, cmdPrint},  //非  (patrolling) symbol for print
    {CMD_CHAR, cmdChar},
    {0x9053, cmdChar},  //道   (wrong) symbol for character
    {CMD_READ, cmdRead},
    {0x2EA2, cmdReead},  //⺢  (way) symbol for read
    {CMD_MOVE, cmdMove},
    {0x5B57, cmdMove,  //字  (water) symbol for move
    {CMD_NULL, cmdNull},
    {0x2AAD, cmdNull},  // 読 (roof+child) symbol for null
    {CMD_PATH, cmdPath},
    {0x2F79, cmdPath},  // 艮 (road) symbol for path
    {CMD_STARTSERVER, cmdStartServer}, 
    {0x2F79, cmdStartserver},  //⻔  (net) symbol for start server
    {CMD_STOPSERVER, cmdStopServer},
    {0x2FA4, cmdDivStopserver},  //⼉  (gate) symbol for stop server
    {CMD_SETATTRIBUTE, cmdSetAttribute},
    {0x2F09, cmdSetattribute},  //⼐  (legs) symbol for set attribute
    {CMD_SETPROXY, cmdSetProxy},
    {0x2F10, cmdDivSetproxy},  // ⼝ (open box) symbol for set proxy
    {CMD_BIN2FUZ, cmdBin2Fuz},
    {0x2F1D, cmdBinfuz},  // ⼍ (mouth) symbol for binary to fuzzy ogic
    {CMD_FUZ2BIN, cmdFuz2Bin},
    {0x2F0D, cmdFuz2bin},  // ⼹ (cover) symbol for fuzzy to binary logic
    {CMD_IMPLY, cmdImply},
    {0x2F39, cmdImply},  //⼺  (snout) symbol for imply
    {CMD_NIMPLY, cmdNimply},
    {0x2F47, cmdNimply},  // ⽇ (bristle) symbol for nimply
    {CMD_NULL, cmdNull},
    {0x2F28, cmdNull},  //⼨ (sun) symbol for null
    {CMD_INHIBITXY, cmdInhibitXY},
    {0x2F6B, cmdInhibitXY},  // ⽫ (inch) symbol for inihbit xy
    {CMD_BUFFER, cmdBuffer},
    {0x2F6B, cmdBuffer},  //⽬ (dish) symbol for buffer
    {CMD_INHIBIT, cmdInhibit},
    {0x2F4A, cmdInhibit},  // ⽊ (eye) symbol for inhibit
    {CMD_EQUIVALENCE, cmdEquivalence}, 
    {0x2F65, cmdEquivalence},  //⽥  (tree) symbol for equivalence
    {CMD_IMPLICATION, cmdImplication}, 
    {0x2F71, cmdImplication},  //⽱  (field) symbol for implication
    {CMD_IDENTITY, cmdIdentity}
    {0x2F7E, cmdIdentity},  //⽾  (track) symbol for identity
    (0x2F7A, cmd), //  ⽺   (plow) symbol for  
    (0x2F60, cmd),  // ⽠    (sheep) symmbol for 
    (0x2F72, cmd),  // ⽲    (melon) symbol for 
     (0x, cmd)           ,  //         (grain) symbol for 
};

#ifdef ENABLE_ADVANCED_MATH
void cmdAdvancedMath(double* operands, int numOperands) {
    // Advanced math operations
};

#else
std::unordered_map<CommandType, CommandFunction> commandMap = {
    {CMD_CMP, cmdCmp},
    {CMD_ADD, cmdAdd},
    // Add all original symbols here...
};

#endif
// building a fuzzy logic and advanced math module:

// Define constants used in logic operations
const double E = 2.718281828459045;
const double FUZZY_100 = 100.0;
const double ZERO = 0.0;
const double ONE = 1.0;
const double FUZZY_50 = 50.0;

// Buffers to simulate SIMD-like vector operations
std::vector<double> vector1 = {0.0, 1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0};
std::vector<double> vector2 = {7.0, 6.0, 5.0, 4.0, 3.0, 2.0, 1.0, 0.0};
std::vector<double> result_vector(8, 0.0);

// Convert binary to fuzzy logic (0 or 1 to 0 or 100)
void binaryToFuzzy(double& value) {
    value = (value != 0) ? FUZZY_100 : ZERO;
}

// Convert fuzzy logic to binary (0-100 to 0 or 1)
void fuzzyToBinary(double& value) {
    value = (value >= FUZZY_50) ? ONE : ZERO;
}

// Compare two vectors and store the result in result_vector
void cmpRoutine(const std::vector<double>& vec1, const std::vector<double>& vec2, std::vector<double>& result) {
    for (size_t i = 0; i < vec1.size(); ++i) {
        double diff = vec1[i] - vec2[i];
        result[i] = (diff > 0) ? diff : 0.0;
    }
}

// Vector AND operation (element-wise minimum)
void vectorAndOptimized(const std::vector<double>& vec1, const std::vector<double>& vec2, std::vector<double>& result) {
    for (size_t i = 0; i < vec1.size(); ++i) {
        result[i] = std::min(vec1[i], vec2[i]);
    }
}

// Fuzzy OR operation (element-wise maximum)
void orRoutine(double& a, double b) {
    a = std::max(a, b);
}

// Fuzzy XOR operation (element-wise difference)
void xorRoutine(double& a, double b) {
    a = std::abs(a - b);
}

// Fuzzy NOT operation (inversion)
void notRoutine(double& a) {
    a = FUZZY_100 - a;
}

// Fuzzy NAND operation (NOT of AND)
void nandRoutine(double& a, double b) {
    a = FUZZY_100 - std::min(a, b);
}

// Fuzzy NOR operation (NOT of OR)
void norRoutine(double& a, double b) {
    a = FUZZY_100 - std::max(a, b);
}

// IF routine for fuzzy logic
void ifRoutine(double& result, double condition, double trueResult, double falseResult) {
    result = (condition > 0) ? trueResult : falseResult;
}

// Implication routine for fuzzy logic
void implicationRoutine(double& a, double b) {
    a = std::min(FUZZY_100, FUZZY_100 - a + b);
}

// Not-implication routine for fuzzy logic
void nimplyRoutine(double& a, double b) {
    a = std::max(0.0, b - a);
}

// NULL routine (reset to zero)
void nullRoutine(double& a) {
    a = 0.0;
}

// Inhibition routine (fuzzy inhibition)
void inhibitionRoutine(double& a, double b) {
    a = std::min(a, FUZZY_100 - b);
}

// Buffer routine (no operation needed in C++ context)
void bufferRoutine(double& a) {
    // No operation
}

// Inhibition XY routine (fuzzy inhibition with two variables)
void inhibitionXYRoutine(double& a, double b) {
    a = std::min(FUZZY_100 - a, b);
}

// Equivalence routine (fuzzy equivalence)
void equivalenceRoutine(double& a, double b) {
    double minAB = std::min(a, b);
    double minNotA = FUZZY_100 - a;
    double minNotB = FUZZY_100 - b;
    a = minAB + std::min(minNotA, minNotB);
}

// Identity routine (always returns 100 in fuzzy logic)
void identityRoutine(double& a) {
    a = FUZZY_100;
}

// Command logic handler for different operations
void logicHandler(char command, double& a, double b) {
    switch (command) {
        case '≈':
            cmpRoutine(vector1, vector2, result_vector);
            break;
        case '&':
            vectorAndOptimized(vector1, vector2, result_vector);
            break;
        case '∨':
            orRoutine(a, b);
            break;
        case '⊕':
            xorRoutine(a, b);
            break;
        case '¬':
            notRoutine(a);
            break;
        case '⊼':
            nandRoutine(a, b);
            break;
        case '⊽':
            norRoutine(a, b);
            break;
        case '⇒':
            ifRoutine(a, a, b, 0); // Example values for trueResult and falseResult
            break;
        case '⇔':
            implicationRoutine(a, b);
            break;
        case '?':
            nimplyRoutine(a, b);
            break;
        case '∄':
            nullRoutine(a);
            break;
        case '⊥':
            inhibitionXYRoutine(a, b);
            break;
        case '∤':
            bufferRoutine(a);
            break;
        case '↔':
            inhibitionRoutine(a, b);
            break;
        case '≡':
            equivalenceRoutine(a, b);
            break;
        case '≃':
            identityRoutine(a);
            break;
        default:
            std::cerr << "Unknown logic command." << std::endl;
            break;
    }
}

// end of code for the advancwed math and fuzzy logic module


// Now continue with all function definitions below...

// Command functions implementation
void cmdCmp(double* operands, int numOperands) {
    if (numOperands >= 2) {
        double result = (operands[0] > operands[1]) ? 100 : 0;  // Fuzzy comparison
        std::cout << "Result (CMP): " << result << std::endl;
    }
}

// Add implementations for all other functions similarly...

// Typedef for function pointer for command execution
typedef std::function<void(double*, int)> CommandFunction;

// Utility function to clamp values to the range [0, 100]
double clamp(double value) {
    return std::max(0.0, std::min(100.0, value));
}

// Main function to execute commands
void executeCommand(CommandType commandType, double* operands, int numOperands) {
    auto it = commandMap.find(commandType);
    if (it != commandMap.end()) {
        it->second(operands, numOperands);
    } else {
        std::cerr << "Unknown command" << std::endl;
        std::cerr << "Error code: " << ERR_UNRECOGNIZED << std::endl;
    }
}

int main() {
    double operands[] = {5.0, 3.0};
    int numOperands = 2;

    executeCommand(CMD_ADD, operands, numOperands);  // Example usage for ADD
    executeCommand(CMD_SUB, operands, numOperands);  // Example usage for SUB

    return 0;
}
