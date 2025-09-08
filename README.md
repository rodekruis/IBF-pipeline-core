# IBF-pipeline-core

This library contains the core functionalities and utilities used across various IBF pipelines. It provides a standardized way to handle common tasks, ensuring consistency and efficiency in pipeline development.

## Main Components

*. `Settings`: Manages configuration settings for the pipelines, from a standardized YAML file. It provides methods to load, validate, and access configuration parameters.
*. `Secrets`: Handles sensitive information such as API keys, passwords, and other credentials. It ensures that secrets are securely stored and accessed only when necessary.
*. `Logger`: Logging utility that provides a standardized way to log messages, errors, and other information during pipeline execution.
*. `Load`: Provides functions to load data from various sources, such as databases or APIs.
*. `Module`: Base class for creating pipeline modules, providing common functionalities and interfaces for module development.




