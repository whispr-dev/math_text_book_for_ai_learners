// `utils.cpp` an attempt to make a C++ version of woflang

#include "utils.h"
#include "err_chk.h" // Include this if you are using error codes in utils.cpp
#include <iostream>
#include <unordered_map>
#include <algorithm> // For std::min and std::max
#include <functional> // For std::function
#include <thread>
#include <mutex>
#include <condition_variable>
#include <atomic>
# include "commands.h"


// Global thread management variables
std::vector<std::thread> threads;
std::mutex mutex_lock;
std::condition_variable cv;
std::atomic<bool> is_thread_ready{false};

// Error handling function
void printError(const char* message) {
    std::cerr << "Error: " << message << std::endl;
}

// Command parsing function with error checking
// CommandType parseCommand(const std::string& commandStr) {
//     static std::unordered_map<std::string, CommandType> commandMap = {
//        {"ADD", CMD_ADD},
//        {"SUB", CMD_SUB},
        // Add more mappings here
//    };

{
    auto it = commandMap.find(commandStr);
    if (it != commandMap.end()) {
        return it->second;
    } else {
        printError("Unknown command"); // Error if command not found
        return CMD_UNKNOWN;  // Return an unknown command type
    }
}

// Math utility function with error checking
double clamp(double value, double min, double max) {
    if (min > max) {
        printError("Invalid clamp range: min is greater than max");
        return value; // No clamping applied if range is invalid
    }
    return std::max(min, std::min(max, value));
}

// Thread management functions with error checking
void createThread(const std::function<void()>& task) {
    try {
        threads.emplace_back(task); // Attempt to create a new thread
    } catch (const std::system_error& e) { // Catch thread creation errors
        printError("Failed to create thread: System error");
    } catch (...) {
        printError("Failed to create thread: Unknown error");
    }
}

void joinAllThreads() {
    for (auto& thread : threads) {
        if (thread.joinable()) {
            try {
                thread.join(); // Attempt to join the thread
            } catch (const std::system_error& e) { // Catch errors when joining threads
                printError("Failed to join thread: System error");
            } catch (...) {
                printError("Failed to join thread: Unknown error");
            }
        } else {
            printError("Thread not joinable or already joined");
        }
    }
}

#include "err_chk.h"

void* safeAllocate(size_t size) {
    void* ptr = malloc(size);
    if (!ptr) {
        reportError(ERR_MEM_INVALID_ADDRESS, "safeAllocate");
        return nullptr;
    }
    return ptr;
}
