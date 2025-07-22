# Exotic Options Pricing - Production Makefile
# ==============================================

# Compiler settings
CXX = g++
CXXFLAGS_BASE = -std=c++17 -Wall -Wextra -Wpedantic -Wconversion -Wshadow
CXXFLAGS_DEBUG = $(CXXFLAGS_BASE) -g -O0 -DDEBUG -DENABLE_MEMORY_PROFILING -fsanitize=address,undefined
CXXFLAGS_RELEASE = $(CXXFLAGS_BASE) -O3 -DNDEBUG -march=native -flto
CXXFLAGS_PROFILE = $(CXXFLAGS_BASE) -O2 -g -pg -DENABLE_PROFILING
CXXFLAGS_TEST = $(CXXFLAGS_BASE) -g -O1 -DENABLE_MEMORY_PROFILING -DENABLE_TESTING

# Directories
SRC_DIR = src
TEST_DIR = tests
BUILD_DIR = build
BIN_DIR = bin
OBJ_DIR = $(BUILD_DIR)/obj
TEST_OBJ_DIR = $(BUILD_DIR)/test_obj

# Source files
SOURCES = $(wildcard $(SRC_DIR)/**/*.cpp) $(wildcard $(SRC_DIR)/*.cpp)
HEADERS = $(wildcard $(SRC_DIR)/**/*.hpp) $(wildcard $(SRC_DIR)/*.hpp)
TEST_SOURCES = $(wildcard $(TEST_DIR)/*.cpp)

# Object files
OBJECTS = $(SOURCES:$(SRC_DIR)/%.cpp=$(OBJ_DIR)/%.o)
TEST_OBJECTS = $(TEST_SOURCES:$(TEST_DIR)/%.cpp=$(TEST_OBJ_DIR)/%.o)

# Targets
MAIN_TARGET = $(BIN_DIR)/exotic_options
TEST_TARGET = $(BIN_DIR)/test_runner

# Libraries
LIBS = -pthread -lm
TEST_LIBS = $(LIBS)

# Create directories
$(shell mkdir -p $(OBJ_DIR)/models $(OBJ_DIR)/utils $(OBJ_DIR)/config)
$(shell mkdir -p $(TEST_OBJ_DIR))
$(shell mkdir -p $(BIN_DIR))

# Default target
.PHONY: all
all: release

# Release build
.PHONY: release
release: CXXFLAGS = $(CXXFLAGS_RELEASE)
release: $(MAIN_TARGET)

# Debug build
.PHONY: debug
debug: CXXFLAGS = $(CXXFLAGS_DEBUG)
debug: LIBS += -fsanitize=address,undefined
debug: $(MAIN_TARGET)

# Profile build
.PHONY: profile
profile: CXXFLAGS = $(CXXFLAGS_PROFILE)
profile: $(MAIN_TARGET)

# Test build
.PHONY: test
test: CXXFLAGS = $(CXXFLAGS_TEST)
test: $(TEST_TARGET)
	@echo "Running unit tests..."
	@./$(TEST_TARGET)

# Build main executable
$(MAIN_TARGET): $(OBJECTS) $(SRC_DIR)/main.cpp
	@echo "Linking $(MAIN_TARGET)..."
	@$(CXX) $(CXXFLAGS) -o $@ $(SRC_DIR)/main.cpp $(OBJECTS) $(LIBS)
	@echo "Build complete: $(MAIN_TARGET)"

# Build test executable
$(TEST_TARGET): $(OBJECTS) $(TEST_OBJECTS)
	@echo "Linking $(TEST_TARGET)..."
	@$(CXX) $(CXXFLAGS) -o $@ $(TEST_OBJECTS) $(filter-out %/main.o,$(OBJECTS)) $(TEST_LIBS)
	@echo "Test build complete: $(TEST_TARGET)"

# Compile source files
$(OBJ_DIR)/%.o: $(SRC_DIR)/%.cpp
	@echo "Compiling $<..."
	@$(CXX) $(CXXFLAGS) -c $< -o $@

# Compile test files
$(TEST_OBJ_DIR)/%.o: $(TEST_DIR)/%.cpp
	@echo "Compiling test $<..."
	@$(CXX) $(CXXFLAGS) -c $< -o $@

# Static analysis
.PHONY: analyze
analyze:
	@echo "Running static analysis..."
	@cppcheck --enable=all --std=c++17 --suppress=missingIncludeSystem $(SRC_DIR)/ $(TEST_DIR)/
	@echo "Static analysis complete"

# Code formatting
.PHONY: format
format:
	@echo "Formatting code..."
	@find $(SRC_DIR) $(TEST_DIR) -name "*.cpp" -o -name "*.hpp" | xargs clang-format -i
	@echo "Code formatting complete"

# Memory leak detection
.PHONY: memcheck
memcheck: debug
	@echo "Running memory leak detection..."
	@valgrind --tool=memcheck --leak-check=full --show-leak-kinds=all --track-origins=yes ./$(MAIN_TARGET)

# Performance profiling
.PHONY: perf
perf: profile
	@echo "Running performance profiling..."
	@./$(MAIN_TARGET)
	@gprof ./$(MAIN_TARGET) gmon.out > performance_report.txt
	@echo "Performance report generated: performance_report.txt"

# Thread safety analysis
.PHONY: thread-check
thread-check: debug
	@echo "Running thread safety analysis..."
	@valgrind --tool=helgrind ./$(MAIN_TARGET)

# Coverage analysis
.PHONY: coverage
coverage: CXXFLAGS = $(CXXFLAGS_TEST) --coverage
coverage: LIBS += --coverage
coverage: clean test
	@echo "Generating coverage report..."
	@gcov $(SOURCES)
	@lcov --capture --directory . --output-file coverage.info
	@genhtml coverage.info --output-directory coverage_html
	@echo "Coverage report generated in coverage_html/"

# Benchmark tests
.PHONY: benchmark
benchmark: release
	@echo "Running benchmarks..."
	@./$(MAIN_TARGET) --benchmark
	@echo "Benchmark complete"

# Integration tests
.PHONY: integration-test
integration-test: release
	@echo "Running integration tests..."
	@python3 -m pytest integration_tests/ -v
	@echo "Integration tests complete"

# Documentation generation
.PHONY: docs
docs:
	@echo "Generating documentation..."
	@doxygen Doxyfile
	@echo "Documentation generated in docs/"

# Install dependencies
.PHONY: install-deps
install-deps:
	@echo "Installing dependencies..."
	@sudo apt-get update
	@sudo apt-get install -y build-essential g++ cmake valgrind cppcheck clang-format lcov doxygen
	@echo "Dependencies installed"

# Package for distribution
.PHONY: package
package: release test
	@echo "Creating distribution package..."
	@mkdir -p dist/exotic_options
	@cp $(MAIN_TARGET) dist/exotic_options/
	@cp README.md dist/exotic_options/
	@cp config.json dist/exotic_options/
	@tar -czf dist/exotic_options.tar.gz -C dist exotic_options
	@echo "Package created: dist/exotic_options.tar.gz"

# Docker build
.PHONY: docker
docker:
	@echo "Building Docker image..."
	@docker build -t exotic-options .
	@echo "Docker image built: exotic-options"

# Clean build artifacts
.PHONY: clean
clean:
	@echo "Cleaning build artifacts..."
	@rm -rf $(BUILD_DIR) $(BIN_DIR)
	@rm -f *.gcov *.gcda *.gcno coverage.info gmon.out
	@rm -rf coverage_html
	@rm -f performance_report.txt
	@rm -f *.log
	@echo "Clean complete"

# Deep clean (including dependencies)
.PHONY: distclean
distclean: clean
	@echo "Deep cleaning..."
	@rm -rf dist/
	@rm -rf docs/
	@echo "Deep clean complete"

# Help target
.PHONY: help
help:
	@echo "Exotic Options Pricing - Build System"
	@echo "====================================="
	@echo ""
	@echo "Available targets:"
	@echo "  all              - Build release version (default)"
	@echo "  release          - Build optimized release version"
	@echo "  debug            - Build debug version with sanitizers"
	@echo "  profile          - Build profiling version"
	@echo "  test             - Build and run unit tests"
	@echo "  analyze          - Run static code analysis"
	@echo "  format           - Format code with clang-format"
	@echo "  memcheck         - Run memory leak detection"
	@echo "  perf             - Run performance profiling"
	@echo "  thread-check     - Run thread safety analysis"
	@echo "  coverage         - Generate code coverage report"
	@echo "  benchmark        - Run performance benchmarks"
	@echo "  integration-test - Run integration tests"
	@echo "  docs             - Generate documentation"
	@echo "  package          - Create distribution package"
	@echo "  docker           - Build Docker image"
	@echo "  install-deps     - Install system dependencies"
	@echo "  clean            - Clean build artifacts"
	@echo "  distclean        - Deep clean including docs/dist"
	@echo "  help             - Show this help message"
	@echo ""
	@echo "Examples:"
	@echo "  make release     - Build optimized version"
	@echo "  make test        - Run all unit tests"
	@echo "  make memcheck    - Check for memory leaks"

# Dependency tracking
-include $(OBJECTS:.o=.d)
-include $(TEST_OBJECTS:.o=.d)

# Generate dependency files
$(OBJ_DIR)/%.d: $(SRC_DIR)/%.cpp
	@$(CXX) $(CXXFLAGS) -MM -MT $(@:.d=.o) $< > $@

$(TEST_OBJ_DIR)/%.d: $(TEST_DIR)/%.cpp
	@$(CXX) $(CXXFLAGS) -MM -MT $(@:.d=.o) $< > $@