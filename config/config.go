// Package config contains configuration constants.
package config

//go:generate go run github.com/fluhus/goat -nh -i config.got -o config.go -df ../config.json

const (
	// DataDir is the main data directory.
	DataDir = ""

	// WSDataDir is the workspace data directory.
	WSDataDir = ""
)
