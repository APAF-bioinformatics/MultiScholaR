# Author(s): Pablo Galaviz
# Email: cmri-bioinformatics@cmri.org.au
# Children’s Medical Research Institute, finding cures for childhood genetic diseases


PROJECT := ProteomeRiver

include config.mk

# ---------------------------------------------------------------
# Directories
# ---------------------------------------------------------------
BASE = $(shell /bin/pwd)
SCRIPTS = $(BASE)/Source
TEST = $(BASE)/tests
SRCD = $(BASE)/R


ifeq ($(UNAME),Darwin)
BASH_PROFILE=$(HOME)/.profile
else
BASH_PROFILE=$(HOME)/.bashrc
endif

intro:
	@echo ''
	@echo '===> $(PROJECT)'
	@echo ''

setdirs:
	@mkdir -p $(PREFIX)/$(PROJECT)/Source/R
	@mkdir -p $(PREFIX)/$(PROJECT)/Source/Python
	@mkdir -p $(PREFIX)/$(PROJECT)/Source/shell
	@mkdir -p $(PREFIX)/$(PROJECT)/test
	@mkdir -p $(PREFIX)/$(PROJECT)/man
	@mkdir -p $(RLIB)


done:
	@echo ''
	@echo '============ All Done! ============'
	@echo ' '


# ---------------------------------------------------------------
# target : help
# ---------------------------------------------------------------

help: intro
	@echo '	'
	@echo '	---------------------------'
	@echo '	Makefile - HELP'
	@echo '	---------------------------'
	@echo '	'
	@echo '	make install			: install $(PROJECT)'
	@echo '	'
	@echo '	make setup			: set environment variable PATH '
	@echo '	'
	@echo '	make uninstall			: uninstall $(PROJECT)'
	@echo ' '


# ---------------------------------------------------------------
# target : install
# ---------------------------------------------------------------

install:intro setdirs
	@echo '	'
	@echo '-----------------------------------'
	@echo 'Installing $(PROJECT) on $(PREFIX)'
	@echo '-----------------------------------'
	@echo '	'
	@$(RSCRIPTEXEC) $(SCRIPTS)/aux/install.R
	@$(REXEC) CMD INSTALL -l $(RLIB) ./
	@cp -r $(SCRIPTS)/R/* $(PREFIX)/$(PROJECT)/Source/R
	@chmod +x $(PREFIX)/$(PROJECT)/Source/R/*
	@cp -r $(SCRIPTS)/Python/* $(PREFIX)/$(PROJECT)/Source/Python
	@chmod +x $(PREFIX)/$(PROJECT)/Source/Python/*
	@cp -r $(TEST)/* $(PREFIX)/$(PROJECT)/test/
	@echo '	'
	@echo 'append $(PREFIX)/$(PROJECT)/Source/R to your PATH environment variable'
	@echo "export PATH=$(PREFIX)/$(PROJECT)/Source/R:$(PREFIX)/$(PROJECT)/Source/Python:\$$PATH"
	@echo "export R_LIBS_USER=$(RLIB)"
	@echo 'or run: make setup'
	@echo '	'
	@echo " === install done ! === "

# ---------------------------------------------------------------
# target : uninstall
# ---------------------------------------------------------------

uninstall:intro
	@echo '	'
	@echo '-----------------------------------'
	@echo 'Uninstalling $(PROJECT) from $(PREFIX)'
	@echo '-----------------------------------'
	@echo '	'
	@rm -rf $(PREFIX)/$(PROJECT)
	@echo '	'
	@echo " === uninstall done ! === "


# ---------------------------------------------------------------
# target : setup
# ---------------------------------------------------------------

setup:intro
	@echo "export PATH=$(PREFIX)/$(PROJECT)/Source/R:\$$PATH" >> $(BASH_PROFILE)
	@echo "export R_LIBS_USER=$(RLIB)" >> $(BASH_PROFILE)
	@echo '	'
	@echo '-----------------------------------'
	@echo "run: source $(BASH_PROFILE)"
	@echo '-----------------------------------'
	@echo '	'
	@echo " === setup done ! === "
