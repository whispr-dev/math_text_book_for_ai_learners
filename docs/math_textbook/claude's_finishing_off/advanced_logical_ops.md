#include <iostream>
#include <vector>
#include <string>
#include <map>
#include <set>
#include <algorithm>
#include <bitset>
#include <sstream>

namespace AdvancedLogicalOperations {
    
    // Represents a boolean expression node
    struct BoolNode {
        enum Type { VAR, NOT, AND, OR, XOR, IMPLIES, EQUIV };
        Type type;
        std::string var;  // for VAR type
        std::vector<BoolNode*> children;
        
        BoolNode(Type t) : type(t) {}
        BoolNode(const std::string& v) : type(VAR), var(v) {}
        
        ~BoolNode() {
            for (auto child : children) delete child;
        }
    };
    
    // Boolean Expression Simplifier using algebraic laws
    class BooleanSimplifier {
    private:
        // Apply De Morgan's laws
        BoolNode* applyDeMorgan(BoolNode* node) {
            if (node->type == BoolNode::NOT && node->children[0]->type == BoolNode::AND) {
                // NOT(A AND B) = NOT(A) OR NOT(B)
                BoolNode* orNode = new BoolNode(BoolNode::OR);
                for (auto child : node->children[0]->children) {
                    BoolNode* notChild = new BoolNode(BoolNode::NOT);
                    notChild->children.push_back(cloneNode(child));
                    orNode->children.push_back(notChild);
                }
                delete node;
                return simplify(orNode);
            }
            if (node->type == BoolNode::NOT && node->children[0]->type == BoolNode::OR) {
                // NOT(A OR B) = NOT(A) AND NOT(B)
                BoolNode* andNode = new BoolNode(BoolNode::AND);
                for (auto child : node->children[0]->children) {
                    BoolNode* notChild = new BoolNode(BoolNode::NOT);
                    notChild->children.push_back(cloneNode(child));
                    andNode->children.push_back(notChild);
                }
                delete node;
                return simplify(andNode);
            }
            return node;
        }
        
        // Remove double negation
        BoolNode* removeDoubleNegation(BoolNode* node) {
            if (node->type == BoolNode::NOT && 
                node->children[0]->type == BoolNode::NOT) {
                BoolNode* result = cloneNode(node->children[0]->children[0]);
                delete node;
                return simplify(result);
            }
            return node;
        }
        
        // Apply absorption laws
        BoolNode* applyAbsorption(BoolNode* node) {
            if (node->type == BoolNode::OR) {
                // A OR (A AND B) = A
                for (size_t i = 0; i < node->children.size(); i++) {
                    for (size_t j = 0; j < node->children.size(); j++) {
                        if (i != j && node->children[j]->type == BoolNode::AND) {
                            if (containsNode(node->children[j], node->children[i])) {
                                // Remove the AND term
                                node->children.erase(node->children.begin() + j);
                                return simplify(node);
                            }
                        }
                    }
                }
            }
            if (node->type == BoolNode::AND) {
                // A AND (A OR B) = A
                for (size_t i = 0; i < node->children.size(); i++) {
                    for (size_t j = 0; j < node->children.size(); j++) {
                        if (i != j && node->children[j]->type == BoolNode::OR) {
                            if (containsNode(node->children[j], node->children[i])) {
                                // Remove the OR term
                                node->children.erase(node->children.begin() + j);
                                return simplify(node);
                            }
                        }
                    }
                }
            }
            return node;
        }
        
        BoolNode* cloneNode(BoolNode* node) {
            if (!node) return nullptr;
            BoolNode* newNode = (node->type == BoolNode::VAR) ? 
                new BoolNode(node->var) : new BoolNode(node->type);
            for (auto child : node->children) {
                newNode->children.push_back(cloneNode(child));
            }
            return newNode;
        }
        
        bool containsNode(BoolNode* haystack, BoolNode* needle) {
            if (haystack->type == needle->type) {
                if (haystack->type == BoolNode::VAR) {
                    return haystack->var == needle->var;
                }
            }
            for (auto child : haystack->children) {
                if (containsNode(child, needle)) return true;
            }
            return false;
        }
        
    public:
        BoolNode* simplify(BoolNode* node) {
            if (!node) return nullptr;
            
            // Recursively simplify children
            for (auto& child : node->children) {
                child = simplify(child);
            }
            
            // Apply simplification rules
            node = removeDoubleNegation(node);
            node = applyDeMorgan(node);
            node = applyAbsorption(node);
            
            // Remove redundant terms in AND/OR
            if (node->type == BoolNode::AND || node->type == BoolNode::OR) {
                // Remove duplicates
                std::set<std::string> seen;
                auto it = node->children.begin();
                while (it != node->children.end()) {
                    if ((*it)->type == BoolNode::VAR) {
                        if (seen.count((*it)->var)) {
                            delete *it;
                            it = node->children.erase(it);
                        } else {
                            seen.insert((*it)->var);
                            ++it;
                        }
                    } else {
                        ++it;
                    }
                }
                
                // If only one child remains, return it
                if (node->children.size() == 1) {
                    BoolNode* child = node->children[0];
                    node->children.clear();
                    delete node;
                    return child;
                }
            }
            
            return node;
        }
        
        std::string toString(BoolNode* node) {
            if (!node) return "";
            
            switch (node->type) {
                case BoolNode::VAR:
                    return node->var;
                case BoolNode::NOT:
                    return "¬" + toString(node->children[0]);
                case BoolNode::AND: {
                    std::string result = "(";
                    for (size_t i = 0; i < node->children.size(); i++) {
                        if (i > 0) result += " ∧ ";
                        result += toString(node->children[i]);
                    }
                    return result + ")";
                }
                case BoolNode::OR: {
                    std::string result = "(";
                    for (size_t i = 0; i < node->children.size(); i++) {
                        if (i > 0) result += " ∨ ";
                        result += toString(node->children[i]);
                    }
                    return result + ")";
                }
                case BoolNode::XOR:
                    return "(" + toString(node->children[0]) + " ⊕ " + 
                           toString(node->children[1]) + ")";
                case BoolNode::IMPLIES:
                    return "(" + toString(node->children[0]) + " → " + 
                           toString(node->children[1]) + ")";
                case BoolNode::EQUIV:
                    return "(" + toString(node->children[0]) + " ↔ " + 
                           toString(node->children[1]) + ")";
            }
            return "";
        }
    };
    
    // Truth Table Generator
    class TruthTableGenerator {
    private:
        std::vector<std::string> variables;
        
        bool evaluate(BoolNode* node, const std::map<std::string, bool>& values) {
            switch (node->type) {
                case BoolNode::VAR:
                    return values.at(node->var);
                case BoolNode::NOT:
                    return !evaluate(node->children[0], values);
                case BoolNode::AND: {
                    bool result = true;
                    for (auto child : node->children) {
                        result = result && evaluate(child, values);
                    }
                    return result;
                }
                case BoolNode::OR: {
                    bool result = false;
                    for (auto child : node->children) {
                        result = result || evaluate(child, values);
                    }
                    return result;
                }
                case BoolNode::XOR:
                    return evaluate(node->children[0], values) != 
                           evaluate(node->children[1], values);
                case BoolNode::IMPLIES:
                    return !evaluate(node->children[0], values) || 
                           evaluate(node->children[1], values);
                case BoolNode::EQUIV:
                    return evaluate(node->children[0], values) == 
                           evaluate(node->children[1], values);
            }
            return false;
        }
        
        void extractVariables(BoolNode* node, std::set<std::string>& vars) {
            if (node->type == BoolNode::VAR) {
                vars.insert(node->var);
            }
            for (auto child : node->children) {
                extractVariables(child, vars);
            }
        }
        
    public:
        std::vector<std::vector<bool>> generate(BoolNode* expression) {
            std::set<std::string> varSet;
            extractVariables(expression, varSet);
            variables = std::vector<std::string>(varSet.begin(), varSet.end());
            
            int n = variables.size();
            int rows = 1 << n;  // 2^n rows
            std::vector<std::vector<bool>> table;
            
            for (int i = 0; i < rows; i++) {
                std::map<std::string, bool> values;
                std::vector<bool> row;
                
                // Set variable values based on binary representation of i
                for (int j = 0; j < n; j++) {
                    bool val = (i >> (n - j - 1)) & 1;
                    values[variables[j]] = val;
                    row.push_back(val);
                }
                
                // Evaluate expression
                bool result = evaluate(expression, values);
                row.push_back(result);
                table.push_back(row);
            }
            
            return table;
        }
        
        void printTable(const std::vector<std::vector<bool>>& table) {
            // Print header
            for (const auto& var : variables) {
                std::cout << var << "\t";
            }
            std::cout << "Result" << std::endl;
            
            // Print separator
            for (size_t i = 0; i <= variables.size(); i++) {
                std::cout << "---\t";
            }
            std::cout << std::endl;
            
            // Print rows
            for (const auto& row : table) {
                for (bool val : row) {
                    std::cout << (val ? "T" : "F") << "\t";
                }
                std::cout << std::endl;
            }
        }
    };
    
    // Karnaugh Map for minimization
    class KarnaughMap {
    private:
        int numVars;
        std::vector<std::vector<int>> map;  // -1: don't care, 0: false, 1: true
        std::vector<std::string> variables;
        
        // Gray code generation
        std::vector<int> grayCode(int n) {
            std::vector<int> result;
            for (int i = 0; i < (1 << n); i++) {
                result.push_back(i ^ (i >> 1));
            }
            return result;
        }
        
        // Check if a rectangle is valid (all 1s or don't cares)
        bool isValidRectangle(int r1, int c1, int r2, int c2) {
            for (int i = r1; i <= r2; i++) {
                for (int j = c1; j <= c2; j++) {
                    if (map[i % map.size()][j % map[0].size()] == 0) {
                        return false;
                    }
                }
            }
            return true;
        }
        
    public:
        KarnaughMap(int vars) : numVars(vars) {
            int rows = 1 << ((vars + 1) / 2);
            int cols = 1 << (vars / 2);
            map = std::vector<std::vector<int>>(rows, std::vector<int>(cols, 0));
        }
        
        void setTruthTable(const std::vector<bool>& truthTable) {
            auto grayRows = grayCode((numVars + 1) / 2);
            auto grayCols = grayCode(numVars / 2);
            
            for (size_t i = 0; i < grayRows.size(); i++) {
                for (size_t j = 0; j < grayCols.size(); j++) {
                    int index = (grayRows[i] << (numVars / 2)) | grayCols[j];
                    if (index < truthTable.size()) {
                        map[i][j] = truthTable[index] ? 1 : 0;
                    }
                }
            }
        }
        
        // Find prime implicants (simplified - basic implementation)
        std::vector<std::string> minimize() {
            std::vector<std::string> terms;
            int rows = map.size();
            int cols = map[0].size();
            std::vector<std::vector<bool>> covered(rows, std::vector<bool>(cols, false));
            
            // Try to find largest rectangles (powers of 2)
            for (int size = rows * cols; size >= 1; size /= 2) {
                for (int r = 0; r < rows; r++) {
                    for (int c = 0; c < cols; c++) {
                        if (map[r][c] == 1 && !covered[r][c]) {
                            // Try different rectangle sizes
                            for (int h = rows; h >= 1; h /= 2) {
                                for (int w = cols; w >= 1; w /= 2) {
                                    if (h * w == size && isValidRectangle(r, c, r + h - 1, c + w - 1)) {
                                        // Mark as covered
                                        for (int i = r; i < r + h; i++) {
                                            for (int j = c; j < c + w; j++) {
                                                covered[i % rows][j % cols] = true;
                                            }
                                        }
                                        // Generate term (simplified representation)
                                        terms.push_back("Group[" + std::to_string(r) + "," + 
                                                      std::to_string(c) + "," + 
                                                      std::to_string(h) + "x" + 
                                                      std::to_string(w) + "]");
                                        goto next_cell;
                                    }
                                }
                            }
                        }
                        next_cell:;
                    }
                }
            }
            
            return terms;
        }
        
        void printMap() {
            std::cout << "Karnaugh Map:" << std::endl;
            for (const auto& row : map) {
                for (int val : row) {
                    if (val == -1) std::cout << "X ";
                    else std::cout << val << " ";
                }
                std::cout << std::endl;
            }
        }
    };
    
    // General demonstration function
    void calc() {
        std::cout << "=== Boolean Algebra Simplification ===" << std::endl;
        
        // Create expression: (A AND B) OR (A AND NOT B)
        BoolNode* expr = new BoolNode(BoolNode::OR);
        BoolNode* and1 = new BoolNode(BoolNode::AND);
        and1->children.push_back(new BoolNode("A"));
        and1->children.push_back(new BoolNode("B"));
        BoolNode* and2 = new BoolNode(BoolNode::AND);
        and2->children.push_back(new BoolNode("A"));
        BoolNode* notB = new BoolNode(BoolNode::NOT);
        notB->children.push_back(new BoolNode("B"));
        and2->children.push_back(notB);
        expr->children.push_back(and1);
        expr->children.push_back(and2);
        
        BooleanSimplifier simplifier;
        std::cout << "Original: " << simplifier.toString(expr) << std::endl;
        BoolNode* simplified = simplifier.simplify(expr);
        std::cout << "Simplified: " << simplifier.toString(simplified) << std::endl;
        
        std::cout << "\n=== Truth Table Generation ===" << std::endl;
        TruthTableGenerator ttGen;
        auto truthTable = ttGen.generate(simplified);
        ttGen.printTable(truthTable);
        
        std::cout << "\n=== Karnaugh Map ===" << std::endl;
        KarnaughMap kmap(2);
        std::vector<bool> ttValues;
        for (const auto& row : truthTable) {
            ttValues.push_back(row.back());
        }
        kmap.setTruthTable(ttValues);
        kmap.printMap();
        auto minimized = kmap.minimize();
        std::cout << "Minimized terms: ";
        for (const auto& term : minimized) {
            std::cout << term << " ";
        }
        std::cout << std::endl;
        
        delete simplified;
    }
}

int main() {
    AdvancedLogicalOperations::calc();
    return 0;
}